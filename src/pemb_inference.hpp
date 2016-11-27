#ifndef PEMB_INFERENCE_HPP
#define PEMB_INFERENCE_HPP

class my_infer {
public:

	static void inference_step(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, gsl_rng *semilla) {
		double logp;

		// 0. Create zero-valued vectors/matrices of counts (for negative samples)
		Matrix1D<int> neg_item_counts = Matrix1D<int>(data.Nitems);
		set_to_zero(neg_item_counts);
		//    Also create auxiliary matrix of negative gradients for alpha
		Matrix2D<double> grad_neg_alpha_aux;
		if(param.flag_userVec==3) {
			grad_neg_alpha_aux = Matrix2D<double>(data.Nitems,param.K);
			set_to_zero(grad_neg_alpha_aux);
		}

		// 1. Compute the sum ot the alpha's
		if(!param.flag_ppca) {
			std::cout << "  computing sum of alphas..." << endl;
			compute_sum_alpha(data,hyper,param,pvar);
		}

		// 2. Initialize all gradients to prior
		std::cout << "  initializing gradients..." << endl;
		if(param.flag_regularization) {
			logp = set_grad_to_prior(data,hyper,param,pvar);
		} else {
			logp = set_grad_to_zero(data,hyper,param,pvar);
		}

		// 3. Increase the gradients with the datapoints where y>0
		std::cout << "  processing positive observations..." << endl;
		logp += process_pos_observations(data,hyper,param,pvar);
		
		// 4. Increase the gradients with (a subset of) the datapoints where y=0
		std::cout << "  processing negative observations..." << endl;
		logp += process_neg_observations(semilla,data,hyper,param,pvar,neg_item_counts,grad_neg_alpha_aux);
		
		// 5. Scale all the negative gradients
		std::cout << "  scaling negative gradients..." << endl;
		scale_neg_gradient(data,hyper,param,pvar,neg_item_counts,grad_neg_alpha_aux);

		// 6. Take gradient step
		std::cout << "  taking grad step..." << endl;
		take_grad_step(data,hyper,param,pvar);

		// 7. Output the objective function to a file
		my_output::write_objective_function(param,logp);
	}

	static double set_grad_to_prior(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		double logp = 0.0;
		logp += set_grad_to_prior_mat(pvar.rho,hyper.s2rho);
		logp += set_grad_to_prior_mat(pvar.alpha,hyper.s2alpha);
		if(param.flag_itemIntercept) {
			logp += set_grad_to_prior_mat(pvar.lambda0,hyper.s2alpha);
		}
		if(param.flag_userVec>0) {
			logp += set_grad_to_prior_mat(pvar.theta,hyper.s2theta);
		}
		if(param.flag_price>0) {
			logp += set_grad_to_prior_mat(pvar.gamma,hyper.s2gamma);
			logp += set_grad_to_prior_mat(pvar.beta,hyper.s2beta);
		}
		return logp;
	}

	static double set_grad_to_prior_mat(Matrix1D<my_pvar_aux> &M, double val) {
		double logp = 0.0;
		for(int i=0; i<M.get_size1(); i++) {
			logp += M.get_object(i).set_grad_to_prior(val);
		}
		return logp;
	}

	static double set_grad_to_prior_mat(Matrix2D<my_pvar_aux> &M, double val) {
		double logp = 0.0;
		for(int i=0; i<M.get_size1(); i++) {
			for(int j=0; j<M.get_size2(); j++) {
				logp += M.get_object(i,j).set_grad_to_prior(val);
			}
		}
		return logp;
	}

	static double set_grad_to_zero(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		set_grad_to_zero_mat(pvar.rho);
		set_grad_to_zero_mat(pvar.alpha);
		if(param.flag_itemIntercept) {
			set_grad_to_zero_mat(pvar.lambda0);
		}
		if(param.flag_userVec>0) {
			set_grad_to_zero_mat(pvar.theta);
		}
		if(param.flag_price>0) {
			set_grad_to_zero_mat(pvar.gamma);
			set_grad_to_zero_mat(pvar.beta);
		}
		return 0.0;
	}

	static void set_grad_to_zero_mat(Matrix1D<my_pvar_aux> &M) {
		for(int i=0; i<M.get_size1(); i++) {
			M.get_object(i).set_grad_to_zero();
		}
	}

	static void set_grad_to_zero_mat(Matrix2D<my_pvar_aux> &M) {
		for(int i=0; i<M.get_size1(); i++) {
			for(int j=0; j<M.get_size2(); j++) {
				M.get_object(i,j).set_grad_to_zero();
			}
		}
	}

	static void compute_sum_alpha(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		// Loop over all transactions
		for(int t=0; t<data.Ntrans; t++) {
			double aux;
			int i;
			double y;
			for(int k=0; k<param.K; k++) {
				aux = 0.0;
				// Loop over all lines of transaction t
				for(const int &ll : data.lines_per_trans.get_object(t)) {
					i = data.obs.y_item[ll];
					if(param.flag_binarizeContext) {
						y = 1.0;
					} else {
						y = static_cast<double>(data.obs.y_rating[ll]);
					}
					aux += y*pvar.alpha.get_object(i,k).e_x;
				}
				pvar.sum_alpha.set_object(t,k,aux);
			}
		}
	}

	static double process_pos_observations(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		double logp;
		int u;
		int s;
		int i;
		int t;
		double y;
		double *suma = new double[param.K];
		double mm;
		double argument;
		double price;
		double dinner_dz;
		double y_minus_mean;
		int j;
		int ll;
		double y_j;
		double denomAvgContext;

		// Initialize logp
		if(param.flag_gaussian) {
			logp = -0.5*static_cast<double>(data.obs.T)*my_log(2.0*M_PI*hyper.s2noise);
		} else if(param.flag_bernoulli) {
			logp = 0.0;
		} else {
			logp = -data.sum_log_yfact;
		}

		// For each line
		for(unsigned int idx_t=0; idx_t<data.obs.T; idx_t++) {
			// Get the user, item, session, transaction, and value
			u = data.obs.y_user[idx_t];
			i = data.obs.y_item[idx_t];
			s = data.obs.y_sess[idx_t];
			t = data.obs.y_trans[idx_t];
			if(param.flag_bernoulli) {
				y = 1.0;
			} else {
				y = static_cast<double>(data.obs.y_rating[idx_t]);
			}

			// Remove the contribution of item i from sum(y_j*alpha_j)
			substract_contribution_sum_alpha(t,i,y,data,param,pvar,suma);

			// If avgContext, scale 'suma' accordingly
			if((param.flag_avgContext==1) && (data.items_per_trans.get_object(t).size()>1)) {
				denomAvgContext = data.Nitems_per_trans.get_object(t)-y;
				for(int k=0; k<param.K; k++) {
					suma[k] /= denomAvgContext;
				}
			} else if((param.flag_avgContext==2) && (data.items_per_trans.get_object(t).size()>1)) {
				denomAvgContext = data.items_per_trans.get_object(t).size()-1.0;
				for(int k=0; k<param.K; k++) {
					suma[k] /= denomAvgContext;
				}
			} else {
				denomAvgContext = 1.0;
			}

			// Compute the Poisson mean
			price = data.get_price(i,s,param);
			mm = compute_mean(param,pvar,i,u,suma,argument,price);
			y_minus_mean = y-mm;

			// Increase log_p
			if(param.flag_bernoulli) {
				logp += my_log(mm);
			} else if(param.flag_gaussian) {
				logp -= 0.5*my_pow2(y_minus_mean)/hyper.s2noise;
			} else if(param.flag_additiveMean) {
				logp += y*my_log(mm)-mm;				
			} else {
				logp += y*argument-mm;				
			}

			// Increase gradient of rho
			for(int k=0; k<param.K; k++) {
				if(param.flag_userVec==0) {
					dinner_dz = suma[k];
				} else if(param.flag_userVec==1) {
					dinner_dz = suma[k]+pvar.theta.get_object(u,k).e_x;
				} else if(param.flag_userVec==2) {
					dinner_dz = suma[k]*pvar.theta.get_object(u,k).e_x;
				} else if(param.flag_userVec==3) {
					dinner_dz = suma[k];
				}
				if(param.flag_gaussian) {
					pvar.rho.get_object(i,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
				} else if(param.flag_bernoulli) {
					pvar.rho.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
				} else if(param.flag_additiveMean) {
					pvar.rho.get_object(i,k).increase_grad(dinner_dz*(y/mm-1.0));
				} else {
					pvar.rho.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
				}
			}
			// Increase gradient of alpha for all items in the context
			for(unsigned int idx_j=0; idx_j<data.items_per_trans.get_object(t).size(); idx_j++) {
				// Find item j (in context)
				j = static_cast<int>(data.items_per_trans.get_object(t).at(idx_j));
				// If item j is indeed in the context (i.e., it is not item i again)
				if(j!=i) {
					// Find y_j (number of units purchased of item j)
					if(param.flag_binarizeContext || param.flag_bernoulli) {
						y_j = 1.0;
					} else {
						ll = data.lines_per_trans.get_object(t).at(idx_j);
						y_j = static_cast<double>(data.obs.y_rating[ll]);
					}
					// Scale down y_j to account for the averaging of elements in context
					y_j /= denomAvgContext;
					// Increase gradient
					for(int k=0; k<param.K; k++) {
						if(param.flag_userVec==0) {
							dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
						} else if(param.flag_userVec==1) {
							dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
						} else if(param.flag_userVec==2) {
							dinner_dz = y_j*pvar.theta.get_object(u,k).e_x*pvar.rho.get_object(i,k).e_x;
						} else if(param.flag_userVec==3) {
							dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
						}
						if(param.flag_gaussian) {
							pvar.alpha.get_object(j,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
						} else if(param.flag_bernoulli) {
							pvar.alpha.get_object(j,k).increase_grad(y_minus_mean*dinner_dz);
						} else if(param.flag_additiveMean) {
							pvar.alpha.get_object(j,k).increase_grad(dinner_dz*(y/mm-1.0));
						} else {
							pvar.alpha.get_object(j,k).increase_grad(y_minus_mean*dinner_dz);
						}
					}
				}
			}
			// Increase gradient of alpha_i (only if param.flag_userVec==3)
			if(param.flag_userVec==3) {
				for(int k=0; k<param.K; k++) {
					dinner_dz = pvar.theta.get_object(u,k).e_x;
					if(param.flag_gaussian) {
						pvar.alpha.get_object(i,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
					} else if(param.flag_bernoulli) {
						pvar.alpha.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
					} else if(param.flag_additiveMean) {
						pvar.alpha.get_object(i,k).increase_grad(dinner_dz*(y/mm-1.0));
					} else {
						pvar.alpha.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
					}
				}
			}
			// Increase gradient of lambda0
			if(param.flag_itemIntercept) {
				if(param.flag_gaussian) {
					pvar.lambda0.get_object(i).increase_grad(y_minus_mean/hyper.s2noise);
				} else if(param.flag_bernoulli) {
					pvar.lambda0.get_object(i).increase_grad(y_minus_mean);
				} else if(param.flag_additiveMean) {
					pvar.lambda0.get_object(i).increase_grad(y/mm-1.0);
				} else {
					pvar.lambda0.get_object(i).increase_grad(y_minus_mean);
				}
			}
			// Increase gradient of theta
			if(param.flag_userVec>0) {
				for(int k=0; k<param.K; k++) {
					if(param.flag_userVec==1) {
						dinner_dz = pvar.rho.get_object(i,k).e_x;
					} else if(param.flag_userVec==2) {
						dinner_dz = suma[k]*pvar.rho.get_object(i,k).e_x;
					} else if(param.flag_userVec==3) {
						dinner_dz = pvar.alpha.get_object(i,k).e_x;
					}
					if(param.flag_gaussian) {
						pvar.theta.get_object(u,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
					} else if(param.flag_bernoulli) {
						pvar.theta.get_object(u,k).increase_grad(y_minus_mean*dinner_dz);
					} else if(param.flag_additiveMean) {
						pvar.theta.get_object(u,k).increase_grad(dinner_dz*(y/mm-1.0));
					} else {
						pvar.theta.get_object(u,k).increase_grad(y_minus_mean*dinner_dz);
					}
				}
			}
			// Increase gradient of the price vectors
			for(int k=0; k<param.flag_price; k++) {
				// gamma
				dinner_dz = -price*pvar.beta.get_object(i,k).e_x;
				if(param.flag_gaussian) {
					pvar.gamma.get_object(u,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
				} else if(param.flag_bernoulli) {
					pvar.gamma.get_object(u,k).increase_grad(y_minus_mean*dinner_dz);
				} else if(param.flag_additiveMean) {
					std::cerr << "[ERR] Additive means cannot be used in combination with price" << endl;
					assert(0);
				} else {
					pvar.gamma.get_object(u,k).increase_grad(y_minus_mean*dinner_dz);
				}
				// beta
				dinner_dz = -price*pvar.gamma.get_object(u,k).e_x;
				if(param.flag_gaussian) {
					pvar.beta.get_object(i,k).increase_grad(y_minus_mean*dinner_dz/hyper.s2noise);
				} else if(param.flag_bernoulli) {
					pvar.beta.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
				} else if(param.flag_additiveMean) {
					std::cerr << "[ERR] Additive means cannot be used in combination with price" << endl;
					assert(0);
				} else {
					pvar.beta.get_object(i,k).increase_grad(y_minus_mean*dinner_dz);
				}
			}
		}
		delete [] suma;
		return logp;
	}

	static double process_neg_observations(gsl_rng *semilla, my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, \
										   Matrix1D<int> &neg_item_counts, Matrix2D<double> &grad_neg_alpha_aux) {
		double logp;
		int u;
		int s;
		int i;
		double *suma = new double[param.K];
		double mm;
		double argument;
		double price;
		double dinner_dz;
		int j;
		int ll;
		double y_j;
		double y_minus_mean;
		double denomAvgContext;
		int count_i;

		// Initialize logp
		if(param.flag_gaussian) {
			logp = -0.5*static_cast<double>(data.Ntrans*param.negsamples)*my_log(2.0*M_PI*hyper.s2noise);
		} else {
			logp = 0.0;
		}

		// For each transaction
		for(int t=0; t<data.Ntrans; t++) {
			// Get the user and session
			u = data.user_per_trans.get_object(t);
			s = data.session_per_trans.get_object(t);

			// Sample negative items from the corresponding distribution
			count_i = 0;
			while(count_i<param.negsamples) {
				i = gsl_ran_discrete(semilla,data.negsampling_dis);
				// If this is indeed a negative sample
				if(data.get_items_in_trans(i,t)<0) {
					// Increase the counts of processed negative samples
					count_i++;
					neg_item_counts.set_object(i,neg_item_counts.get_object(i)+1);
					// Set 'suma' to sum_alpha
					set_sum_alpha(t,data,param,pvar,suma);
					// If avgContext, scale 'suma' accordingly
					if(param.flag_avgContext==1) {
						denomAvgContext = data.Nitems_per_trans.get_object(t);
						for(int k=0; k<param.K; k++) {
							suma[k] /= denomAvgContext;
						}
					} else if(param.flag_avgContext==2) {
						denomAvgContext = data.items_per_trans.get_object(t).size();
						for(int k=0; k<param.K; k++) {
							suma[k] /= denomAvgContext;
						}
					} else {
						denomAvgContext = 1.0;
					}
					// Compute the Poisson mean
					price = data.get_price(i,s,param);
					mm = compute_mean(param,pvar,i,u,suma,argument,price);
					y_minus_mean = -mm;
					// Increase log_p
					if(param.flag_bernoulli) {
						logp += my_log(1.0-mm);
					} else if(param.flag_gaussian) {
						logp -= 0.5*my_pow2(mm)/hyper.s2noise;
					} else if(param.flag_additiveMean) {
						logp -= mm;
					} else {
						logp -= mm;
					}
					// Increase gradient of rho
					for(int k=0; k<param.K; k++) {
						if(param.flag_userVec==0) {
							dinner_dz = suma[k];
						} else if(param.flag_userVec==1) {
							dinner_dz = suma[k]+pvar.theta.get_object(u,k).e_x;
						} else if(param.flag_userVec==2) {
							dinner_dz = suma[k]*pvar.theta.get_object(u,k).e_x;
						} else if(param.flag_userVec==3) {
							dinner_dz = suma[k];
						}
						if(param.flag_gaussian) {
							pvar.rho.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz/hyper.s2noise);
						} else if(param.flag_bernoulli) {
							pvar.rho.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz);
						} else if(param.flag_additiveMean) {
							pvar.rho.get_object(i,k).increase_grad_neg(-dinner_dz);
						} else {
							pvar.rho.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz);
						}
					}
					// Increase gradient of alpha for all items in the context
					for(unsigned int idx_j=0; idx_j<data.items_per_trans.get_object(t).size(); idx_j++) {
						// Find item j (in context)
						j = static_cast<int>(data.items_per_trans.get_object(t).at(idx_j));
						// Find y_j (number of units purchased of item j)
						if(param.flag_binarizeContext || param.flag_bernoulli) {
							y_j = 1.0;
						} else {
							ll = data.lines_per_trans.get_object(t).at(idx_j);
							y_j = static_cast<double>(data.obs.y_rating[ll]);
						}
						// Scale down y_j to account for the averaging of elements in context
						y_j /= denomAvgContext;
						// Increase gradient
						for(int k=0; k<param.K; k++) {
							if(param.flag_userVec==0) {
								dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
							} else if(param.flag_userVec==1) {
								dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
							} else if(param.flag_userVec==2) {
								dinner_dz = y_j*pvar.theta.get_object(u,k).e_x*pvar.rho.get_object(i,k).e_x;
							} else if(param.flag_userVec==3) {
								dinner_dz = y_j*pvar.rho.get_object(i,k).e_x;
							}
							if(param.flag_gaussian) {
								pvar.alpha.get_object(j,k).increase_grad_neg(y_minus_mean*dinner_dz/hyper.s2noise);
							} else if(param.flag_bernoulli) {
								pvar.alpha.get_object(j,k).increase_grad_neg(y_minus_mean*dinner_dz);
							} else if(param.flag_additiveMean) {
								pvar.alpha.get_object(j,k).increase_grad_neg(-dinner_dz);
							} else {
								pvar.alpha.get_object(j,k).increase_grad_neg(y_minus_mean*dinner_dz);
							}
						}
					}
					// Increase gradient of alpha_i (only if param.flag_userVec==3)
					if(param.flag_userVec==3) {
						double aux;
						for(int k=0; k<param.K; k++) {
							dinner_dz = pvar.theta.get_object(u,k).e_x;
							aux = grad_neg_alpha_aux.get_object(i,k);
							if(param.flag_gaussian) {
								grad_neg_alpha_aux.set_object(i,k,aux+y_minus_mean*dinner_dz/hyper.s2noise);
							} else if(param.flag_bernoulli) {
								grad_neg_alpha_aux.set_object(i,k,aux+y_minus_mean*dinner_dz);
							} else if(param.flag_additiveMean) {
								grad_neg_alpha_aux.set_object(i,k,aux-dinner_dz);
							} else {
								grad_neg_alpha_aux.set_object(i,k,aux+y_minus_mean*dinner_dz);
							}
						}
					}
					// Increase gradient of lambda0
					if(param.flag_itemIntercept) {
						if(param.flag_gaussian) {
							pvar.lambda0.get_object(i).increase_grad_neg(y_minus_mean/hyper.s2noise);
						} else if(param.flag_bernoulli) {
							pvar.lambda0.get_object(i).increase_grad_neg(y_minus_mean);
						} else if(param.flag_additiveMean) {
							pvar.lambda0.get_object(i).increase_grad_neg(-1.0);
						} else {
							pvar.lambda0.get_object(i).increase_grad_neg(y_minus_mean);
						}
					}
					// Increase gradient of theta
					if(param.flag_userVec>0) {
						for(int k=0; k<param.K; k++) {
							if(param.flag_userVec==1) {
								dinner_dz = pvar.rho.get_object(i,k).e_x;
							} else if(param.flag_userVec==2) {
								dinner_dz = suma[k]*pvar.rho.get_object(i,k).e_x;
							} else if(param.flag_userVec==3) {
								dinner_dz = pvar.alpha.get_object(i,k).e_x;
							}
							if(param.flag_gaussian) {
								pvar.theta.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz/hyper.s2noise);
							} else if(param.flag_bernoulli) {
								pvar.theta.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz);
							} else if(param.flag_additiveMean) {
								pvar.theta.get_object(u,k).increase_grad_neg(-dinner_dz);
							} else {
								pvar.theta.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz);
							}
						}
					}
					// Increase gradient of the price vectors
					for(int k=0; k<param.flag_price; k++) {
						// gamma
						dinner_dz = -price*pvar.beta.get_object(i,k).e_x;
						if(param.flag_gaussian) {
							pvar.gamma.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz/hyper.s2noise);
						} else if(param.flag_bernoulli) {
							pvar.gamma.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz);
						} else if(param.flag_additiveMean) {
							std::cerr << "[ERR] Additive means cannot be used in combination with price" << endl;
							assert(0);
						} else {
							pvar.gamma.get_object(u,k).increase_grad_neg(y_minus_mean*dinner_dz);
						}
						// beta
						dinner_dz = -price*pvar.gamma.get_object(u,k).e_x;
						if(param.flag_gaussian) {
							pvar.beta.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz/hyper.s2noise);
						} else if(param.flag_bernoulli) {
							pvar.beta.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz);
						} else if(param.flag_additiveMean) {
							std::cerr << "[ERR] Additive means cannot be used in combination with price" << endl;
							assert(0);
						} else {
							pvar.beta.get_object(i,k).increase_grad_neg(y_minus_mean*dinner_dz);
						}
					}
				}
			}
		}
		delete [] suma;
		// Scale up and down logp
		logp *= param.zeroFactor*static_cast<double>(data.Ntrans*data.Nitems-data.obs.T)/static_cast<double>(data.Ntrans*param.negsamples);
		return logp;
	}

	static double compute_mean(const my_param &param, my_pvar &pvar, int i, int u, double *suma, double &argument, double price) {
		double inner_prod = 0.0;
		// Compute inner product rho*suma
		for(int k=0; k<param.K; k++) {
			if(param.flag_userVec==0) {
				inner_prod += pvar.rho.get_object(i,k).e_x*suma[k];
			} else if(param.flag_userVec==1) {
				inner_prod += pvar.rho.get_object(i,k).e_x*(suma[k]+pvar.theta.get_object(u,k).e_x);
			} else if(param.flag_userVec==2) {
				inner_prod += pvar.rho.get_object(i,k).e_x*suma[k]*pvar.theta.get_object(u,k).e_x;
			} else if(param.flag_userVec==3) {
				inner_prod += pvar.rho.get_object(i,k).e_x*suma[k]+pvar.theta.get_object(u,k).e_x*pvar.alpha.get_object(i,k).e_x;
			} else {
				std::cerr << "[ERR] Wrong value of userVec" << endl;
				assert(0);
			}
		}
		// Add intercept term
		if(param.flag_itemIntercept) {
			inner_prod += pvar.lambda0.get_object(i).e_x;
		}
		for(int k=0; k<param.flag_price; k++) {
			inner_prod -= pvar.gamma.get_object(u,k).e_x*pvar.beta.get_object(i,k).e_x*price;
		}
		// Extra computations (e.g., truncate or exponentiate)
		if(param.flag_additiveMean) {
			inner_prod = ((inner_prod<1e-6)?1e-6:inner_prod);
			argument = inner_prod;
		} else if(param.flag_bernoulli) {
			inner_prod = ((inner_prod<-23.0259)?-23.0259:inner_prod);
			inner_prod = ((inner_prod>23.0259)?23.0259:inner_prod);
			argument = inner_prod;
			inner_prod = my_sigmoid(inner_prod);
		} else if(!param.flag_gaussian) {
			argument = inner_prod;
			inner_prod = my_exp(inner_prod);
		}
		// Return
		return inner_prod;
	}

	static void substract_contribution_sum_alpha(int t, int i, double y, my_data &data, const my_param &param, my_pvar &pvar, double *vec) {
		if(param.flag_binarizeContext) {
			y = 1.0;
		}
		// Remove contribution from sum_lambdas
		for(int k=0; k<param.K; k++) {
			vec[k] = pvar.sum_alpha.get_object(t,k)-y*pvar.alpha.get_object(i,k).e_x;
		}
	}

	static void set_sum_alpha(int t, my_data &data, const my_param &param, my_pvar &pvar, double *vec) {
		for(int k=0; k<param.K; k++) {
			vec[k] = pvar.sum_alpha.get_object(t,k);
		}
	}

	static void scale_neg_gradient(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, \
								   Matrix1D<int> &neg_item_counts, Matrix2D<double> &grad_neg_alpha_aux) {
		double sviFactor;
		// per-user variables
		if(param.flag_userVec>0 || param.flag_price>0) {
			for(int u=0; u<data.Nusers; u++) {
				// The sviFactor is: ((#trips of user u)*(#items)-(total #items purchased by user u))/((#trips of user u)*(#negsamples))
				sviFactor = param.zeroFactor*static_cast<double>(data.sessions_per_user.get_object(u).size()*data.Nitems-data.lines_per_user.get_object(u).size())/static_cast<double>(data.sessions_per_user.get_object(u).size()*param.negsamples);
				if(param.flag_userVec>0) {
					for(int k=0; k<param.K; k++) {
						pvar.theta.get_object(u,k).scale_add_grad_neg(sviFactor);
					}
				}
				for(int k=0; k<param.flag_price; k++) {
					pvar.gamma.get_object(u,k).scale_add_grad_neg(sviFactor);
				}
			}
		}
		// per-item variables (except 'alpha')
		for(int i=0; i<data.Nitems; i++) {
			if(neg_item_counts.get_object(i)>0) {
				// The sviFactor is: (#trans-(total #times item i is purchased))/(#negsamples for item i)
				sviFactor = param.zeroFactor*static_cast<double>(data.Ntrans-data.lines_per_item.get_object(i).size())/static_cast<double>(neg_item_counts.get_object(i));
				for(int k=0; k<param.K; k++) {
					pvar.rho.get_object(i,k).scale_add_grad_neg(sviFactor);
				}
				for(int k=0; k<param.flag_price; k++) {
					pvar.beta.get_object(i,k).scale_add_grad_neg(sviFactor);
				}
				if(param.flag_itemIntercept) {
					pvar.lambda0.get_object(i).scale_add_grad_neg(sviFactor);
				}
				if(param.flag_userVec==3) {
					for(int k=0; k<param.K; k++) {
						pvar.alpha.get_object(i,k).increase_grad(sviFactor*grad_neg_alpha_aux.get_object(i,k));
					}
				}
			}
		}
		// alpha
		for(int i=0; i<data.Nitems; i++) {
			// The sviFactor is: ((#trans in which item i appears)*(#items)-(sum_{t: i appears} size(t)))/((#trans in which item i appears)*(#param.negsamples))
			sviFactor = param.zeroFactor*static_cast<double>(data.trans_per_item.get_object(i).size()*data.Nitems-data.sum_sizetrans_per_item.get_object(i))/static_cast<double>(data.trans_per_item.get_object(i).size()*param.negsamples);
			for(int k=0; k<param.K; k++) {
				pvar.alpha.get_object(i,k).scale_add_grad_neg(sviFactor);
			}
		}
	}

	static void compute_test_performance(int duration, my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		if(param.noTest) {
			return;
		}

		double llh = 0.0;
		double llh_fair = 0.0;
		double llh_fair_binary = 0.0;
		int count = 0;
		int count_nitems = 0;
		int u;
		int i;
		int s;
		double y;
		double poiss_mean;
		double *suma = new double[param.K];
		double *p_item = new double[data.Nitems];
		double sum_norm;
		double argument = 0;
		double price = 0.0;

		// Compute sum_alpha for each transaction in test set
		Matrix2D<double> test_sum_alpha = Matrix2D<double>(data.test_Ntrans,param.K);
		set_to_zero(test_sum_alpha);
		for(int t=0; t<data.test_Ntrans; t++) {
			double aux;
			// For each item in transaction t
			for(int &ll : data.test_lines_per_trans.get_object(t)) {
				i = data.obs_test.y_item[ll];
				if(param.flag_binarizeContext) {
					y = 1.0;
				} else {
					y = static_cast<double>(data.obs_test.y_rating[ll]);
				}

				// Only add contribution of items found in train.tsv
				if(data.lines_per_item.get_object(i).size()>0) {
					for(int k=0; k<param.K; k++) {
						aux = test_sum_alpha.get_object(t,k);
						aux += y*pvar.alpha.get_object(i,k).e_x;
						test_sum_alpha.set_object(t,k,aux);
					}
				}
			}
		}

		// Create a vector specifying how many valid lines per transactions
		// there are in test.tsv. A line is valid if both the item and the user
		// appear in train.tsv (the user is necessary only when -userVec != 0)
		Matrix1D<int> test_valid_lines_per_trans = Matrix1D<int>(data.test_Ntrans);
		set_to_zero(test_valid_lines_per_trans);
		for(unsigned int ll=0; ll<data.obs_test.T; ll++) {
			u = data.obs_test.y_user[ll];
			i = data.obs_test.y_item[ll];
			int t = data.obs_test.y_trans[ll];
			if(data.lines_per_item.get_object(i).size()>0) {
				if(param.flag_userVec==0 || data.lines_per_user.get_object(u).size()>0) {
					int aux = 1+test_valid_lines_per_trans.get_object(t);
					test_valid_lines_per_trans.set_object(t,aux);
				}
			}
		}

		// Compute inner product for each line in test set
		count = 0;
		for(unsigned int ll=0; ll<data.obs_test.T; ll++) {
			u = data.obs_test.y_user[ll];
			i = data.obs_test.y_item[ll];
			s = data.obs_test.y_sess[ll];
			int t = data.obs_test.y_trans[ll];
			y = static_cast<double>(data.obs_test.y_rating[ll]);
			double y_rem = y;
			if(param.flag_binarizeContext) {
				y_rem = 1.0;
			}
			price = data.get_price(i,s,param);

			// Ignore "non-valid" transactions
			if(test_valid_lines_per_trans.get_object(t)>1) {
				// Remove the contribution of item i from sum_alpha
				for(int k=0; k<param.K; k++) {
					suma[k] = test_sum_alpha.get_object(t,k)-y_rem*pvar.alpha.get_object(i,k).e_x;
					// If avgContext, scale 'suma' accordingly
					if(param.flag_avgContext==1) {
						suma[k] /= (data.test_Nitems_per_trans.get_object(t)-y_rem);
					}
					else if(param.flag_avgContext==2) {
						suma[k] /= (data.test_items_per_trans.get_object(t).size()-1.0);
					}
				}
				poiss_mean = compute_mean(param,pvar,i,u,suma,argument,price);
				// (A) Compute the likelihood according to the model
				if(param.flag_gaussian) {
					llh -= 0.5*my_log(2.0*M_PI*hyper.s2noise)+0.5*my_pow2(y-poiss_mean)/hyper.s2noise;
				} else if(param.flag_bernoulli) {
					llh += y*my_log(poiss_mean)+(1.0-y)*my_log(1-poiss_mean);
				} else if(param.flag_additiveMean) {
					llh += y*my_log(poiss_mean)-poiss_mean-my_logfactorial(y);
				} else {
					llh += y*argument-poiss_mean-my_logfactorial(y);
				}
				// (B) Compute the "fair" llh (requires iterate over all items j)
				sum_norm = 0.0;
				for(int j=0; j<data.Nitems; j++) {
					price = data.get_price(j,s,param);
					p_item[j] = compute_mean(param,pvar,j,u,suma,argument,price);
					sum_norm += p_item[j];
				}
				// Normalize p_item
				for(int j=0; j<data.Nitems; j++) {
					p_item[j] /= sum_norm;
				}
				// Compute llh_fair
				if(!param.flag_gaussian) {
					llh_fair += y*my_log(p_item[i]);
					llh_fair_binary += my_log(p_item[i]);
				}
				// Increase count to compute average llh later
				count++;
				count_nitems += data.obs_test.y_rating[ll];
			}
		}

		// Take the average
		llh /= static_cast<double>(count);
		llh_fair /= static_cast<double>(count_nitems);
		llh_fair_binary /= static_cast<double>(count);

		// Print to file
		string fname = param.outdir+"/test.txt";
		char buffer[500];
    	sprintf(buffer,"%d\t%d\t%.9f\t%.9f\t%.9f\t%d\t%d",param.it,duration,llh,llh_fair,llh_fair_binary,count,count_nitems);
        my_output::write_line(fname,string(buffer));

		// Free memory
		delete [] suma;
		delete [] p_item;
	}

	static double compute_val_likelihood(bool writeToFile, int duration, my_data &data, const my_param &param, const my_hyper &hyper, my_pvar &pvar) {
		if(param.noVal) {
			return 0.0;
		}

		double llh = 0.0;
		double llh_fair = 0.0;
		double llh_fair_binary = 0.0;
		int count = 0;
		int count_nitems = 0;
		int u;
		int i;
		int s;
		int t;
		double y;
		double poiss_mean;
		double *suma = new double[param.K];
		double *p_item = new double[data.Nitems];
		double sum_norm;
		double argument = 0;
		double price = 0.0;

		// Ensure that sum_alpha is up-to-date
		compute_sum_alpha(data,hyper,param,pvar);

		// For each line of validation.tsv
		for(unsigned int ll=0; ll<data.obs_val.T; ll++) {
			u = data.obs_val.y_user[ll];
			i = data.obs_val.y_item[ll];
			s = data.obs_val.y_sess[ll];
			t = data.obs_val.y_trans[ll];
			y = static_cast<double>(data.obs_val.y_rating[ll]);

			// If the transaction and the item are found in train.tsv
			if(t>=0 && data.lines_per_item.get_object(i).size()>0) {
				// Set 'suma' to sum_alpha
				set_sum_alpha(t,data,param,pvar,suma);
				// If avgContext, scale 'suma' accordingly
				if(param.flag_avgContext==1) {
					for(int k=0; k<param.K; k++) {
						suma[k] /= static_cast<double>(data.Nitems_per_trans.get_object(t));
					}
				} else if(param.flag_avgContext==2) {
					for(int k=0; k<param.K; k++) {
						suma[k] /= static_cast<double>(data.items_per_trans.get_object(t).size());
					}
				}
				// Compute the Poisson mean
				price = data.get_price(i,s,param);
				poiss_mean = compute_mean(param,pvar,i,u,suma,argument,price);
				// (A) Compute the likelihood according to the model
				if(param.flag_gaussian) {
					llh -= 0.5*my_log(2.0*M_PI*hyper.s2noise)+0.5*my_pow2(y-poiss_mean)/hyper.s2noise;
				} else if(param.flag_bernoulli) {
					llh += y*my_log(poiss_mean)+(1.0-y)*my_log(1.0-poiss_mean);
				} else if(param.flag_additiveMean) {
					poiss_mean = (poiss_mean<1e-10)?1e-10:poiss_mean;
					llh += y*my_log(poiss_mean)-poiss_mean-my_logfactorial(y);
				} else {
					poiss_mean = (poiss_mean<1e-10)?1e-10:poiss_mean;
					llh += y*argument-poiss_mean-my_logfactorial(y);
				}
				// (B) Compute the "fair" llh (requires iterate over all items j)
				sum_norm = 0.0;
				for(int j=0; j<data.Nitems; j++) {
					price = data.get_price(j,s,param);
					p_item[j] = compute_mean(param,pvar,j,u,suma,argument,price);
					sum_norm += p_item[j];
				}
				// Normalize p_item
				for(int j=0; j<data.Nitems; j++) {
					p_item[j] /= sum_norm;
				}
				// Compute llh_fair
				if(!param.flag_gaussian) {
					llh_fair += y*my_log(p_item[i]);
					llh_fair_binary += my_log(p_item[i]);
				}
				// Increase count to compute average llh later
				count++;
				count_nitems += data.obs_val.y_rating[ll];
			}
		}

		// Take the average
		llh /= static_cast<double>(count);
		llh_fair /= static_cast<double>(count_nitems);
		llh_fair_binary /= static_cast<double>(count);

		// Print to file
		if(writeToFile) {
			string fname = param.outdir+"/validation.txt";
			char buffer[500];
	    	sprintf(buffer,"%d\t%d\t%.9f\t%.9f\t%.9f\t%d\t%d",param.it,duration,llh,llh_fair,llh_fair_binary,count,count_nitems);
	        my_output::write_line(fname,string(buffer));
		}

		// Free memory
		delete [] suma;
		delete [] p_item;

		// Return
		return llh_fair;
	}

	static double compute_avg_norm(Matrix2D<my_pvar_aux> &m) {
		double suma = 0.0;
		for(int i=0; i<m.get_size1(); i++) {
			for(int j=0; j<m.get_size2(); j++) {
				suma += my_pow2(m.get_object(i,j).e_x);
			}
		}
		return(suma/static_cast<double>(m.get_size1()*m.get_size2()));
	}

	static void set_to_zero(Matrix1D<int> &M) {
		for(int i=0; i<M.get_size1(); i++) {
			M.set_object(i,0);
		}
	}

	static void set_to_zero(Matrix2D<double> &M) {
		for(int i=0; i<M.get_size1(); i++) {
			for(int j=0; j<M.get_size2(); j++) {
				M.set_object(i,j,0.0);
			}
		}
	}

	static void take_grad_step(const my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		for(int i=0; i<data.Nitems; i++) {
			for(int k=0; k<param.K; k++) {
				pvar.rho.get_object(i,k).take_grad_step(param,param.eta,param.gamma,param.it);
				pvar.alpha.get_object(i,k).take_grad_step(param,param.eta,param.gamma,param.it);
			}
			if(param.flag_itemIntercept) {
				pvar.lambda0.get_object(i).take_grad_step(param,param.eta,param.gamma,param.it);
			}
		}
		if(param.flag_userVec>0) {
			for(int u=0; u<data.Nusers; u++) {
				for(int k=0; k<param.K; k++) {
					pvar.theta.get_object(u,k).take_grad_step(param,param.eta,param.gamma,param.it);
				}
			}
		}
		for(int k=0; k<param.flag_price; k++) {
			for(int u=0; u<data.Nusers; u++) {
				pvar.gamma.get_object(u,k).take_grad_step(param,param.eta,param.gamma,param.it);
			}
			for(int i=0; i<data.Nitems; i++) {
				pvar.beta.get_object(i,k).take_grad_step(param,param.eta,param.gamma,param.it);
			}
		}
	}
};

#endif
