#ifndef PEMB_INFERENCE_PPCA_HPP
#define PEMB_INFERENCE_PPCA_HPP

class my_infer_ppca {
public:
	static void inference_step(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, \
							   gsl_rng **semilla_n, int it) {
		// Compute inner products for all ratings
		compute_all_exp_inner_prods(data,hyper,param,pvar);

		// Compute the gradients wrt the embedding vectors
		compute_grad_rho(data,hyper,param,pvar,semilla_n);

		// Compute the gradients wrt the context vectors
		compute_grad_alpha(data,hyper,param,pvar,semilla_n);

		// Take gradient step
		take_grad_step(data,param,pvar,it);
	}

	static void compute_test_performance(int duration, my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		if(param.noTest) {
			return;
		} else {
			std::cerr << "[ERR] Test performance not implemented because" << endl \
					  << "      it involves further optimization" << endl;
			assert(0);
		}
		return;
	}

	static double compute_val_likelihood(bool writeToFile, int duration, my_data &data, const my_param &param, const my_hyper &hyper, my_pvar &pvar) {
		if(param.noVal) {
			return 0;
		}

		double llh = 0.0;
		double llh_fair = 0.0;
		double llh_fair_binary = 0.0;
		int count = 0;
		int count_nitems = 0;
		int i;
		int t;
		double y;
		double poiss_mean;
		double *suma = new double[param.K];
		double *p_item = new double[data.Nitems];
		double sum_norm;

		// For each line of validation.tsv
		for(unsigned int ll=0; ll<data.obs_val.T; ll++) {
			i = data.obs_val.y_item[ll];
			y = static_cast<double>(data.obs_val.y_rating[ll]);

			// Get the transaction (from (user,session) pair)
			t = data.obs_val.y_trans[ll];
			if(t>=0) {
				// Compute the Poisson mean
				poiss_mean = compute_exp_inner_prod(param,pvar,t,i);
				// (A) Compute the likelihood according to the model
				poiss_mean = (poiss_mean<1e-10)?1e-10:poiss_mean;
				llh += y*my_log(poiss_mean)-poiss_mean-my_logfactorial(y);
				// (B) Compute the "fair" llh (requires iterate over all items j)
				sum_norm = 0.0;
				for(int j=0; j<data.Nitems; j++) {
					p_item[j] = compute_exp_inner_prod(param,pvar,t,j);
					sum_norm += p_item[j];
				}
				// Normalize p_item
				for(int j=0; j<data.Nitems; j++) {
					p_item[j] /= sum_norm;
				}
				// Compute llh_fair
				llh_fair += y*my_log(p_item[i]);
				llh_fair_binary += my_log(p_item[i]);
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

	static void compute_grad_rho(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, \
							     gsl_rng **semilla_n) {
		// Initialize with the gradient of the prior
		std::vector<std::thread> threads;
		int thread_incr = data.Ntrans/param.Nthreads;
		int thread_extra = data.Ntrans%param.Nthreads;
		int thread_start = 0;
		int thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int t=t_start; t<t_end; t++) {
					for(int k=0; k<param.K; k++) {
						pvar.rho.get_object(t,k).set_grad_to_prior(hyper.s2rho,param);
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}

		// Main loop
		threads.clear();
		thread_incr = data.Ntrans/param.Nthreads;
		thread_extra = data.Ntrans%param.Nthreads;
		thread_start = 0;
		thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int t=t_start; t<t_end; t++) {
					double y;
					double dif;
					double aux;
					int count;
					int cc;
					int maxNegItems;
					double sviFactor;

					// (A) Loop over items that *appear* in transaction t
					cc = 0;
					for(int &i : data.items_per_trans.get_object(t)) {
						// Find number of units of item i purchased in transaction t
						y = data.get_items_in_trans(i,t);
						// Obtain gradient
						dif = y-pvar.sum_alpha.get_object(t,cc);
						for(int k=0; k<param.K; k++) {
							aux = dif*pvar.alpha.get_object(i,k).e_x;
							pvar.rho.get_object(t,k).increase_grad(aux);
						}
						cc++;
					}
					// (B) Loop over some other items that *do not* appear in transaction t
					unsigned int Nitems_neg = data.Nitems-data.items_per_trans.get_object(t).size();
					// Maximum number of negative samples
					if(param.negsamples<=0) {
						maxNegItems = Nitems_neg;
					} else {
						maxNegItems = my_min(Nitems_neg,param.negsamples);
					}
					// Compute sviFactor
					sviFactor = static_cast<double>(Nitems_neg)/static_cast<double>(maxNegItems);
					if(param.zeroFactor>=0) {
						sviFactor *= param.zeroFactor;
					}
					// Create a vector containing the items over which the "negative" for loop should not iterate
					// Initialize this vector to the list of "positives"
					std::vector<int> items_visited(data.items_per_trans.get_object(t));
					count = 0;
					while(count<maxNegItems) {
						// Generate a random item
						int i = gsl_rng_uniform_int(semilla_n[t_id],data.Nitems);
						// If this item is indeed "negative" _and_ it has not been already visited
						if(!my_infer::is_transaction_visited(items_visited,i)) {
							// Compute exponentiated inner product X sviFactor
							dif = -sviFactor*compute_exp_inner_prod(param,pvar,t,i);
							// Obtain gradient
							for(int k=0; k<param.K; k++) {
								aux = dif*pvar.alpha.get_object(i,k).e_x;
								pvar.rho.get_object(t,k).increase_grad(aux);
							}
							// Increase count of negative samples
							// Update variables for the 'while' loop
							items_visited.push_back(i);
							count++;	
						}
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}
	}

	static void compute_grad_alpha(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar, \
							       gsl_rng **semilla_n) {
		// Initialize with the gradient of the prior
		std::vector<std::thread> threads;
		int thread_incr = data.Nitems/param.Nthreads;
		int thread_extra = data.Nitems%param.Nthreads;
		int thread_start = 0;
		int thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int i=t_start; i<t_end; i++) {
					for(int k=0; k<param.K; k++) {
						pvar.alpha.get_object(i,k).set_grad_to_prior(hyper.s2alpha,param);
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}

		// Main loop
		threads.clear();
		thread_incr = data.Nitems/param.Nthreads;
		thread_extra = data.Nitems%param.Nthreads;
		thread_start = 0;
		thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int i=t_start; i<t_end; i++) {
					double y;
					double dif;
					int t;
					double aux;
					int aux_i;
					int count;
					double sviFactor;
					int maxPosTrans;
					int maxNegTrans;

					// (A) Transactions in which item i *does* appear
					unsigned int Ntrans_pos = data.trans_per_item.get_object(i).size();
					unsigned int *order_pos = new unsigned int[Ntrans_pos];
					for(unsigned int idx_t=0; idx_t<Ntrans_pos; idx_t++) {
						order_pos[idx_t] = idx_t;
					}
					// Shuffle the vector
					gsl_ran_shuffle(semilla_n[t_id],order_pos,Ntrans_pos,sizeof(unsigned int));
					// Maximum number of positive samples
					if(param.possamples<=0) {
						maxPosTrans = Ntrans_pos;
					} else {
						maxPosTrans = my_min(Ntrans_pos,param.possamples);
					}
					// Compute sviFactor
					sviFactor = static_cast<double>(Ntrans_pos)/static_cast<double>(maxPosTrans);
					// Loop
					for(int idx_t=0; idx_t<maxPosTrans; idx_t++) {
						// Loop over (a subset of) the positive transactions
						t = data.trans_per_item.get_object(i).at(order_pos[idx_t]);
						// Find number of units of item i purchased in transaction t
						y = data.get_items_in_trans(i,t);
						// Find the order in which item i appears in transcation t
						aux_i = data.get_order_in_trans(i,t);
						// Obtain gradient
						dif = sviFactor*(y-pvar.sum_alpha.get_object(t,aux_i));
						for(int k=0; k<param.K; k++) {
							aux = dif*pvar.rho.get_object(t,k).e_x;
							pvar.alpha.get_object(i,k).increase_grad(aux);
						}
						// Increase index
						idx_t++;
					}
					// (B) Transactions in which item i *does not* appear
					unsigned int Ntrans_neg = data.Ntrans-Ntrans_pos;
					// Maximum number of negative samples
					if(param.negsamples<=0) {
						maxNegTrans = Ntrans_neg;
					} else {
						maxNegTrans = my_min(Ntrans_neg,param.negsamples);
					}
					// Compute sviFactor
					sviFactor = static_cast<double>(Ntrans_neg)/static_cast<double>(maxNegTrans);
					if(param.zeroFactor>=0) {
						sviFactor *= param.zeroFactor;
					}
					// Create a vector containing the transactions over which the "negative" for loop should not iterate
					// Initialize this vector to the list of "positives"
					std::vector<int> trans_visited(data.trans_per_item.get_object(i));
					count = 0;
					while(count<maxNegTrans) {
						// Generate a random transaction
						t = gsl_rng_uniform_int(semilla_n[t_id],data.Ntrans);
						// If this transaction is indeed "negative" _and_ it has not been visited
						if(!my_infer::is_transaction_visited(trans_visited,t)) {
							// Compute exponentiated inner product X sviFactor
							dif = -sviFactor*compute_exp_inner_prod(param,pvar,t,i);
							// Obtain gradient
							for(int k=0; k<param.K; k++) {
								aux = dif*pvar.rho.get_object(t,k).e_x;
								pvar.alpha.get_object(i,k).increase_grad(aux);
							}
						}
						// Update variables for the 'while' loop
						trans_visited.push_back(t);
						count++;
					}
					// Free memory
					delete [] order_pos;
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}		
	}

	static void compute_all_exp_inner_prods(my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
		std::vector<std::thread> threads;
		int thread_incr = data.Ntrans/param.Nthreads;
		int thread_extra = data.Ntrans%param.Nthreads;
		int thread_start = 0;
		int thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int t=t_start; t<t_end; t++) {
					double aux;
					int cc = 0;
					// For each item in transaction t
					for(int &i : data.items_per_trans.get_object(t)) {
						aux = compute_exp_inner_prod(param,pvar,t,i);
						pvar.sum_alpha.set_object(t,cc,aux);
						cc++;
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}
	}

	inline static double compute_exp_inner_prod(const my_param &param, my_pvar &pvar, int t, int i) {
		double aux = 0.0;
		for(int k=0; k<param.K; k++) {
			aux += pvar.rho.get_object(t,k).e_x*pvar.alpha.get_object(i,k).e_x;
		}
		aux = (aux<-18.4207?-18.4207:aux);
		aux = (aux>3.9120?3.9120:aux);
		return my_exp(aux);
	}

	static void take_grad_step(const my_data &data, const my_param &param, my_pvar &pvar, int it) {
		// Rho
		std::vector<std::thread> threads;
		int thread_incr = data.Ntrans/param.Nthreads;
		int thread_extra = data.Ntrans%param.Nthreads;
		int thread_start = 0;
		int thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int t=t_start; t<t_end; t++) {
					for(int k=0; k<param.K; k++) {
						pvar.rho.get_object(t,k).take_grad_step(param,param.etaRho,param.gammaRho,it);
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}

		// Alpha
		threads.clear();
		thread_incr = data.Nitems/param.Nthreads;
		thread_extra = data.Nitems%param.Nthreads;
		thread_start = 0;
		thread_end = thread_incr;

		for(unsigned int nthread=0; nthread<param.Nthreads; nthread++) {
			if(nthread==param.Nthreads-1) {
				thread_end += thread_extra;
			}
			threads.push_back(std::thread([&](int t_start, int t_end, unsigned int t_id) -> void
			{
				for(int i=t_start; i<t_end; i++) {
					for(int k=0; k<param.K; k++) {
						pvar.alpha.get_object(i,k).take_grad_step(param,param.etaAlpha,param.gammaAlpha,it);
					}
				}
			},thread_start,thread_end,nthread));
			thread_start = thread_end;
			thread_end = thread_start+thread_incr;
		}

		for(std::thread &t_ : threads) {
			t_.join();
		}
	}
};

#endif
