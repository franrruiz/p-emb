#ifndef PEMB_IO_HPP
#define PEMB_IO_HPP

class my_input {
public:

	static void read_input_from_command_line(int argc, char **argv, my_data &data, my_param &param, my_hyper &hyper) {
	    int i = 0;
	    string val;
	    while(i<=argc-1) {
	        if(strcmp(argv[i], "-dir") == 0) {
	            val = string(argv[++i]);
	            val = remove_final_slash(val);
	            param.datadir = val;
	        } else if(strcmp(argv[i], "-outdir") == 0) {
	            val = string(argv[++i]);
	            val = remove_final_slash(val);
	            param.outdir = val;
	        } else if(strcmp(argv[i], "-K") == 0) {
	            val = string(argv[++i]);
	            param.K = std::stoi(val);
	        } else if(strcmp(argv[i], "-seed") == 0) {
	            val = string(argv[++i]);
	        	param.seed = std::stoi(val);
	        } else if(strcmp(argv[i], "-rfreq") == 0) {
	            val = string(argv[++i]);
	        	param.rfreq = std::stoi(val);
	        } else if(strcmp(argv[i], "-saveCycle") == 0) {
	            val = string(argv[++i]);
	        	param.saveCycle = std::stoi(val);
	        } else if(strcmp(argv[i], "-max-iterations") == 0) {
	            val = string(argv[++i]);
	        	param.Niter = std::stoi(val);
	        } else if(strcmp(argv[i], "-label") == 0) {
	            val = string(argv[++i]);
	            param.label = val;
	        } else if(strcmp(argv[i], "-nsFreq") == 0) {
	        	param.flag_nsFreq = true;
	        } else if(strcmp(argv[i], "-negsamples") == 0) {
	            val = string(argv[++i]);
	        	param.negsamples = std::stoi(val);
	        } else if(strcmp(argv[i], "-zeroFactor") == 0) {
	            val = string(argv[++i]);
	        	param.zeroFactor = std::stod(val);
	        } else if(strcmp(argv[i], "-userVec") == 0) {
	            val = string(argv[++i]);
	        	param.flag_userVec = std::stoi(val);
	        } else if(strcmp(argv[i], "-avgContext") == 0) {
	            val = string(argv[++i]);
	        	param.flag_avgContext = std::stoi(val);
	        } else if(strcmp(argv[i], "-itemIntercept") == 0) {
	        	param.flag_itemIntercept = true;
	        } else if(strcmp(argv[i], "-price") == 0) {
	            val = string(argv[++i]);
	        	param.flag_price = std::stoi(val);
	        } else if(strcmp(argv[i], "-normPrice") == 0) {
	        	param.flag_normPrice = true;
	        } else if(strcmp(argv[i], "-additiveMean") == 0) {
	        	param.flag_additiveMean = true;
	        } else if(strcmp(argv[i], "-ppca") == 0) {
	        	param.flag_ppca = true;
	        	std::cout << "[ERR] PPCA is not implemented yet" << endl;
	        	assert(0);
	        } else if(strcmp(argv[i], "-gaussian") == 0) {
	        	param.flag_gaussian = true;
	        } else if(strcmp(argv[i], "-bernoulli") == 0) {
	        	param.flag_bernoulli = true;
	        } else if(strcmp(argv[i], "-rmsprop") == 0) {
	        	param.flag_rmsprop = true;
	        	param.flag_adagrad = false;
	        } else if(strcmp(argv[i], "-adagrad") == 0) {
	        	param.flag_rmsprop = false;
	        	param.flag_adagrad = true;
	        } else if(strcmp(argv[i], "-binarizeContext") == 0) {
	        	param.flag_binarizeContext = true;
	        } else if(strcmp(argv[i], "-noVal") == 0) {
	        	param.noVal = true;
	        } else if(strcmp(argv[i], "-noTest") == 0) {
	        	param.noTest = true;
	        } else if(strcmp(argv[i], "-eta") == 0) {
	            val = string(argv[++i]);
	            param.eta = std::stod(val);
	        } else if(strcmp(argv[i], "-gamma") == 0) {
	            val = string(argv[++i]);
	            param.gamma = std::stod(val);
	        } else if(strcmp(argv[i], "-stdIni") == 0) {
	            val = string(argv[++i]);
	            param.stdIni = std::stod(val);
	        } else if(strcmp(argv[i], "-iniPath") == 0) {
	            val = string(argv[++i]);
	            param.iniPath = remove_final_slash(val);
	        } else if(strcmp(argv[i], "-iniThetaVal") == 0) {
	            val = string(argv[++i]);
	            param.flag_iniThetaVal = true;
	            param.iniThetaVal = std::stod(val);
	        } else if(strcmp(argv[i], "-iniPriceVal") == 0) {
	            val = string(argv[++i]);
	            param.flag_iniPriceVal = true;
	            param.iniPriceVal = std::stod(val);
	        } else if(strcmp(argv[i], "-valTolerance") == 0) {
	            val = string(argv[++i]);
	            param.valTolerance = std::stod(val);
	        } else if(strcmp(argv[i], "-valConsecutive") == 0) {
	            val = string(argv[++i]);
	            param.valConsecutive = std::stoi(val);
	        } else if(strcmp(argv[i], "-noRegulariz") == 0) {
	            param.flag_regularization = false;
	        } else if(strcmp(argv[i], "-s2rho") == 0) {
	            val = string(argv[++i]);
	            hyper.s2rho = std::stod(val);
	        } else if(strcmp(argv[i], "-s2alpha") == 0) {
	            val = string(argv[++i]);
	            hyper.s2alpha = std::stod(val);
	        } else if(strcmp(argv[i], "-s2theta") == 0) {
	            val = string(argv[++i]);
	            hyper.s2theta = std::stod(val);
	        } else if(strcmp(argv[i], "-s2gamma") == 0) {
	            val = string(argv[++i]);
	            hyper.s2gamma = std::stod(val);
	        } else if(strcmp(argv[i], "-s2beta") == 0) {
	            val = string(argv[++i]);
	            hyper.s2beta = std::stod(val);
	        } else if(strcmp(argv[i], "-s2all") == 0) {
	            val = string(argv[++i]);
	            hyper.s2rho = std::stod(val);
	            hyper.s2alpha = std::stod(val);
	            hyper.s2theta = std::stod(val);
	            hyper.s2beta = std::stod(val);
	            hyper.s2gamma = std::stod(val);
	        } else if(strcmp(argv[i], "-s2noise") == 0) {
	            val = string(argv[++i]);
	            hyper.s2noise = std::stod(val);
	        } else if(i>0) {
	            fprintf(stdout,"[ERR] Unknown option %s\n", argv[i]);
	            assert(0);
	        } 
	        i++;
	    };
	    // Check for incompatibilities
	    if(param.flag_gaussian && param.flag_additiveMean) {
	    	std::cerr << "[WARN] Disabling '-additiveMean' because '-gaussian' was found (it is implicit)" << endl;
	    	param.flag_additiveMean = false;
	    }
	    if(param.flag_gaussian && !param.noVal && param.valTolerance>0) {
	    	std::cerr << "[WARN] Setting valTolerance=0 because '-gaussian' was found" << endl;
	    	param.valTolerance = 0.0;
	    }
	    if(param.flag_gaussian && param.flag_bernoulli) {
	    	std::cerr << "[ERR] '-gaussian' and '-bernoulli' cannot be used together" << endl;
	    	assert(0);
	    }
	    if(param.flag_bernoulli && param.flag_additiveMean) {
	    	std::cerr << "[WARN] Disabling '-additiveMean' because '-bernoulli' was found" << endl;
	    	param.flag_additiveMean = false;
	    }
	    if(param.flag_ppca && param.flag_gaussian) {
	    	std::cerr << "[WARN] Disabling '-gaussian' because '-ppca' was found" << endl;
	    	param.flag_gaussian = false;
	    }
	    if(param.flag_ppca && param.flag_bernoulli) {
	    	std::cerr << "[WARN] Disabling '-bernoulli' because '-ppca' was found" << endl;
	    	param.flag_bernoulli = false;
	    }
	    if(param.flag_ppca && param.flag_additiveMean) {
	    	std::cerr << "[WARN] Disabling '-additiveMean' because '-ppca' was found" << endl;
	    	param.flag_additiveMean = false;
	    }
	    if(param.flag_ppca && param.flag_userVec>0) {
	    	std::cerr << "[WARN] Setting userVec=0 because '-ppca' was found" << endl;
	    	param.flag_userVec = 0;
	    }
	    if(param.flag_ppca && !param.noTest) {
	    	std::cerr << "[WARN] Setting noTest=true because '-ppca' was found" << endl;
	    	param.noTest = true;
	    }
	    if(param.flag_ppca && param.flag_itemIntercept) {
	    	std::cerr << "[WARN] Setting itemIntercept=false because '-ppca' was found" << endl;
	    	param.flag_itemIntercept = false;
	    }
	    if(param.flag_ppca && param.flag_price>0) {
	    	std::cerr << "[WARN] Setting price=false because '-ppca' was found" << endl;
	    	param.flag_price = 0;
	    }
	    if(param.flag_binarizeContext && param.flag_avgContext==1) {
	    	std::cerr << "[WARN] Setting avgContext=2 because '-binarizeContext' was found" << endl;
	    	param.flag_avgContext = 2;
	    }
	    if(param.flag_price>0 && param.flag_additiveMean) {
	    	std::cerr << "[WARN] Ignoring '-price' because '-additive' was found" << endl;
	    	param.flag_price = 0;
	    }
	    if(param.flag_nsFreq!=-1 && param.flag_nsFreq!=0 && param.flag_nsFreq!=1) {
	    	std::cerr << "[ERR] Wrong argument for -nsFreq" << endl;
	    	assert(0);
	    }
	    // Set stdIni
        if(param.stdIni<=0.0) {
            param.stdIni = 0.1/sqrt(param.K);
        }
	}

    static string remove_final_slash(const string f) {
        string aux;
        if(f.back()=='/') {
            aux = f.substr(0,f.size()-1);
        } else {
            aux = f;
        }
        return aux;
    }

    static void read_data_file(my_data &data, const my_param &param) {
		string trainFname = param.datadir+"/train.tsv";
		string testFname = param.datadir+"/test.tsv";
		string validationFname = param.datadir+"/validation.tsv";
		unsigned int NiTr;
		unsigned int NuTr;

		// This is a 2-step process
		// Step 1: Count #users, #items (and #sessions)
	    std::cout << " +Stage 1/2..." << endl;
		int n_u = 0;	// Auxiliary index to last user that was read
		int n_i = 0;	// Auxiliary index to last item that was read
		int n_s = 0;	// Auxiliary index to last session that was read
		data.obs.T = count_input_tsv_file(trainFname,data,param,n_u,n_i,n_s);
		NuTr = data.user_ids.size();
		NiTr = data.item_ids.size();
		if(!param.noTest) {
			data.obs_test.T = count_input_tsv_file(testFname,data,param,n_u,n_i,n_s);
		}
		if(!param.noVal) {
			data.obs_val.T = count_input_tsv_file(validationFname,data,param,n_u,n_i,n_s);
		}
		if(data.user_ids.size()>NuTr) {
			std::cerr << "[WARN] There are users in test/validation that do not appear in train" << endl;
		}
		if(data.item_ids.size()>NiTr) {
			std::cerr << "[WARN] There are items in test/validation that do not appear in train" << endl;
		}
		data.Nusers = data.user_ids.size();
		data.Nitems = data.item_ids.size();
		data.Nsessions = data.session_ids.size();

		// Allocate memory to store data
		data.obs = my_data_aux(data.obs.T);
		data.obs_test = my_data_aux(data.obs_test.T);
		data.obs_val = my_data_aux(data.obs_val.T);

		// Step 2: Read the actual data
	    std::cout << " +Stage 2/2..." << endl;
		read_input_tsv_file(trainFname,data,data.obs,param);
		if(!param.noTest) {
			read_input_tsv_file(testFname,data,data.obs_test,param);
		}
		if(!param.noVal) {
			read_input_tsv_file(validationFname,data,data.obs_val,param);
		}
	}

	static void read_a_line(const my_param &param, FILE *fin, unsigned long long *uid, unsigned long long *mid, unsigned long long *sid, unsigned int *rating) {
		fscanf(fin,"%llu\t%llu\t%llu\t%u\n",uid,mid,sid,rating);
	}

	static unsigned int count_input_tsv_file(string fname, my_data &data, const my_param &param, int &n_u, int &n_i, int &n_s) {
		unsigned long long uid;
		unsigned long long mid;
		unsigned long long sid;
		unsigned int rating;
		unsigned int nlines = 0;		// number of lines read

	  	FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
		while(!feof(fin)){
			// Read a line
			read_a_line(param,fin,&uid,&mid,&sid,&rating);
			nlines++;
			// Append user id to the list of user id's
			if(data.user_ids.find(uid) == data.user_ids.end()) {
				data.user_ids.insert(std::pair<unsigned long long,int>(uid,n_u));
				n_u++;
			}
			// Append item id to the list of item id's
			if(data.item_ids.find(mid) == data.item_ids.end()) {
				data.item_ids.insert(std::pair<unsigned long long,int>(mid,n_i));
				n_i++;
			}
			// Append session id to the list of session id's
			if(data.session_ids.find(sid) == data.session_ids.end()) {
				data.session_ids.insert(std::pair<unsigned long long,int>(sid,n_s));
				n_s++;
			}
		}
		fclose(fin);
		return(nlines);
	}

	static void read_input_tsv_file(string fname, const my_data &data, my_data_aux &data_aux, const my_param &param) {
		unsigned long long uid;
		unsigned long long mid;
		unsigned long long sid;
		unsigned int rating;
		unsigned int count = 0;

	  	FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}

		while(!feof(fin)){
			// Read a line
			read_a_line(param,fin,&uid,&mid,&sid,&rating);
			// Store values
			data_aux.y_user[count] = data.user_ids.find(uid)->second;
			data_aux.y_item[count] = data.item_ids.find(mid)->second;
			if(param.flag_bernoulli) {
				data_aux.y_rating[count] = (rating>0)?1:0;
			} else {
				data_aux.y_rating[count] = rating;
			}
			data_aux.y_sess[count] = data.session_ids.find(sid)->second;
			// Increase count
			count++;
		}
		fclose(fin);
	}

    static void read_price_file(my_data &data, const my_param &param) {
    	// Open file
    	string fname = param.datadir+"/item_sess_price.tsv";
		FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
	  	// Allocate data structure
    	data.price_is = Matrix2D<double>(data.Nitems,data.Nsessions);
    	for(int i=0; i<data.Nitems; i++) {
	    	for(int s=0; s<data.Nsessions; s++) {
	    		data.price_is.set_object(i,s,0.0);
	    	}
    	}
    	// Read file
    	int nlines = 0;
    	int ll = 0;
    	unsigned long long mid;
    	unsigned long long sid;
    	double p;
    	std::map<unsigned long long,int>::const_iterator i_item;
    	std::map<unsigned long long,int>::const_iterator i_session;
    	int idx_item;
    	int idx_session;
		while(!feof(fin)){
			ll++;
			fscanf(fin,"%llu\t%llu\t%lf\n",&mid,&sid,&p);
			// Find the item and session
			i_item = data.item_ids.find(mid);
			i_session = data.session_ids.find(sid);
			// Check if the ids are found in the data files
			if(i_item==data.item_ids.end() || i_session==data.session_ids.end()) {
				std::cerr << "[WARN] Line " << ll << " of item_sess_price.tsv contains a session or item id not found in the data files" << endl;
			} else {
				nlines++;
				idx_item = i_item->second;
				idx_session = i_session->second;
				data.price_is.set_object(idx_item,idx_session,p);
			}
		}
    	// Close file
    	fclose(fin);
    	// Check number of valid lines
    	if(nlines<data.Nitems*data.Nsessions) {
    		std::cerr << "[ERR] The number of lines in 'item_sess_price.tsv'";
    		std::cerr << "      is not consistent with the data files" << endl;
    		assert(0);
    	}
    }	
    
    static void normalize_prices(my_data &data, const my_param &param) {
    	data.price_avg_i = Matrix1D<double>(data.Nitems);
    	double aux;

    	for(int i=0; i<data.Nitems; i++) {
    		// Add up all prices of item i
    		aux = 0.0;
	    	for(int s=0; s<data.Nsessions; s++) {
	    		aux += data.price_is.get_object(i,s)/static_cast<double>(data.Nsessions);
	    	}
	    	data.price_avg_i.set_object(i,aux);
    		// Divide prices by their avg
    		if(param.flag_normPrice) {
		    	for(int s=0; s<data.Nsessions; s++) {
		    		aux = data.price_is.get_object(i,s)/data.price_avg_i.get_object(i);
		    		data.price_is.set_object(i,s,aux);
		    	}
    		}
    	}
    }
    
    static void take_log_prices(my_data &data, const my_param &param) {
    	double aux;
    	for(int i=0; i<data.Nitems; i++) {
	    	for(int s=0; s<data.Nsessions; s++) {
	    		aux = data.price_is.get_object(i,s);
	    		data.price_is.set_object(i,s,my_log(aux));
	    	}
	    }
	}
};

class my_output {
public:

	static void create_output_folder(const my_data &data, my_param &param) {
		ostringstream sa;
		sa << "t" << data.Ntrans << "-";
		sa << "n" << data.Nusers << "-";
		sa << "m" << data.Nitems << "-";
		sa << "k" << param.K;
		if(param.flag_itemIntercept) {
			sa << "-intercept";
		}
		sa << "-avgCtxt" << param.flag_avgContext;
		if(param.flag_binarizeContext) {
			sa << "-binCtxt";
		}
		if(param.flag_additiveMean) {
			sa << "-additive";
		}
		if(param.flag_ppca) {
			sa << "-ppca";
		}
		if(param.flag_gaussian) {
			sa << "-gauss";
		}
		if(param.flag_bernoulli) {
			sa << "-bern";
		}
		if(param.flag_userVec>0) {
			sa << "-users" << param.flag_userVec;
		}
		if(param.flag_price>0) {
			sa << "-price" << param.flag_price;
			if(param.flag_normPrice) {
				sa << "norm";
			}
		}
		sa << "-eta" << param.eta;
		if(param.zeroFactor>=0) {
			sa << "-zF" << param.zeroFactor;
		}
		if(param.negsamples>=0) {
			sa << "-nS" << param.negsamples;
		}
		if(param.flag_nsFreq!=-1) {
			sa << "-nsFreq" << param.flag_nsFreq;
		}
		if(param.label!="") {
			sa << "-" << param.label;
		}

		string prefix = sa.str();
		struct stat buffer;
		std::cout << "Output directory label: " << prefix << endl;
		int out = stat((param.outdir+"/"+prefix).c_str(), &buffer);
		if(out != 0) {
			std::cout << "Creating directory " << param.outdir << "/" << prefix << endl;
			if(mkdir((param.outdir+"/"+prefix).c_str(), S_IRUSR | S_IWUSR | S_IXUSR) != 0) {
				std::cerr << "[ERR] Cannot create output folder" << endl;
				assert(0);
			}
		}

		// Keep the new outdir
		param.outdir = param.outdir+"/"+prefix;
	}

    static void create_log_file(const my_data &data, const my_param &param, const my_hyper &hyper) {
        // Write data
        write_log(param.outdir,"Data:");
        write_log(param.outdir," +datadir="+param.datadir);
        write_log(param.outdir," +Nusers="+std::to_string(data.Nusers));
        write_log(param.outdir," +Nitems="+std::to_string(data.Nitems));
        write_log(param.outdir," +Nsessions="+std::to_string(data.Nsessions));
        write_log(param.outdir," +Ntrans="+std::to_string(data.Ntrans));
        write_log(param.outdir," +Ntrans (test)="+std::to_string(data.test_Ntrans));
        write_log(param.outdir," +Lines of train.tsv="+std::to_string(data.obs.T));
        write_log(param.outdir," +Lines of test.tsv="+std::to_string(data.obs_test.T));
        write_log(param.outdir," +Lines of validation.tsv="+std::to_string(data.obs_val.T));

        // Write parameters
        write_log(param.outdir,"Parameters:");
        write_log(param.outdir," +outdir="+param.outdir);
        write_log(param.outdir," +K="+std::to_string(param.K));
        write_log(param.outdir," +seed="+std::to_string(param.seed));
        write_log(param.outdir," +rfreq="+std::to_string(param.rfreq));
        write_log(param.outdir," +saveCycle="+std::to_string(param.saveCycle));
        write_log(param.outdir," +max-iterations="+std::to_string(param.Niter));
        write_log(param.outdir," +negsamples="+std::to_string(param.negsamples));
        write_log(param.outdir," +nsFreq="+std::to_string(param.flag_nsFreq));
        write_log(param.outdir," +zeroFactor="+std::to_string(param.zeroFactor));
        write_log(param.outdir," +userVec="+std::to_string(param.flag_userVec));
        write_log(param.outdir," +avgContext="+std::to_string(param.flag_avgContext));
        write_log(param.outdir," +itemIntercept="+std::to_string(param.flag_itemIntercept));
        write_log(param.outdir," +price="+std::to_string(param.flag_price));
        write_log(param.outdir," +normPrice="+std::to_string(param.flag_normPrice));
        write_log(param.outdir," +binarizeContext="+std::to_string(param.flag_binarizeContext));
        write_log(param.outdir," +additiveMean="+std::to_string(param.flag_additiveMean));
        write_log(param.outdir," +ppca="+std::to_string(param.flag_ppca));
        write_log(param.outdir," +gaussian="+std::to_string(param.flag_gaussian));
        write_log(param.outdir," +bernoulli="+std::to_string(param.flag_bernoulli));
        write_log(param.outdir," +regularization="+std::to_string(param.flag_regularization));
        write_log(param.outdir," +rmsprop="+std::to_string(param.flag_rmsprop));
        write_log(param.outdir," +adagrad="+std::to_string(param.flag_adagrad));
        write_log(param.outdir," +eta="+std::to_string(param.eta));
        write_log(param.outdir," +gamma="+std::to_string(param.gamma));
        write_log(param.outdir," +valTolerance="+std::to_string(param.valTolerance));
        write_log(param.outdir," +valConsecutive="+std::to_string(param.valConsecutive));

        // Write initialization
        write_log(param.outdir,"Initialization:");
        write_log(param.outdir," +stdIni="+std::to_string(param.stdIni));
        write_log(param.outdir," +iniPath="+param.iniPath);
        if(param.flag_iniThetaVal) {
	        write_log(param.outdir," +iniThetaVal="+std::to_string(param.iniThetaVal));
        }
        if(param.flag_iniPriceVal) {
	        write_log(param.outdir," +iniPriceVal="+std::to_string(param.iniPriceVal));
        }

        // Write hyperparameters
        write_log(param.outdir,"Hyperparameters:");
        write_log(param.outdir," +s2rho="+std::to_string(hyper.s2rho));
        write_log(param.outdir," +s2alpha="+std::to_string(hyper.s2alpha));
        write_log(param.outdir," +s2theta="+std::to_string(hyper.s2theta));
        write_log(param.outdir," +s2gamma="+std::to_string(hyper.s2gamma));
        write_log(param.outdir," +s2beta="+std::to_string(hyper.s2beta));
        write_log(param.outdir," +s2noise="+std::to_string(hyper.s2noise));
	}

    static void write_log(const string folder, const string str) {
        write_line(folder+"/log.txt",str);
    }

    static void write_line(const string file, const string str) {
        ofstream myfile;
        myfile.open(file,ios::out|ios::app);    // Append content to the end of existing file
        if(!myfile.is_open()) {
            std::cerr << "[ERR] File '" << file << "' could not be opened" << endl;
            assert(0);
        }
        myfile << str << endl;
        myfile.close();
    }

    static void write_normalized_prices(my_data &data, const my_param &param) {
    	string filename = param.outdir+"/mean_prices.txt";
    	char buffer[200];
    	unsigned long long iName;
    	int i;
    	string aux;

    	int count = 0;
    	for(auto const &it : data.item_ids) {
    		iName = it.first;
    		i = it.second;
    		sprintf(buffer,"%d\t%llu\t%.12f",count,iName,data.price_avg_i.get_object(i));
    		aux = string(buffer);
    		write_line(filename,aux);
    		count++;
    	}
	}

    static void write_telapsed(const string folder, int iter, time_t t_ini, time_t t_end) {
        double elapsed_secs = difftime(t_end,t_ini);
        write_line(folder+"/telapsed.txt",std::to_string(iter)+"\t"+std::to_string(elapsed_secs));
    }

    static void write_norm(const string outdir, int iter, double norm_rho, double norm_alpha) {
        write_line(outdir+"/avg_norm.txt",std::to_string(iter)+"\t"+std::to_string(norm_rho)+"\t"+std::to_string(norm_alpha));
    }

    static void write_norm(const string outdir, int iter, double norm_rho, double norm_alpha, double norm_theta) {
        write_line(outdir+"/avg_norm.txt",std::to_string(iter)+"\t"+std::to_string(norm_rho)+"\t"+std::to_string(norm_alpha)+"\t"+std::to_string(norm_theta));
    }

    static void write_max_file(const my_param &param, int duration, double val_llh, int why) {
    	char buffer[200];
    	sprintf(buffer,"%d\t%d\t%.9f\t%d",param.it,duration,val_llh,why);
        write_line(param.outdir+"/max.txt",string(buffer));
    }
    
    static void write_matrix(string filename, std::map<unsigned long long,int> &ids, Matrix1D<my_pvar_aux> &M) {
    	char buffer[200];
    	string aux;
    	int i;
    	unsigned long long iName;
    	std::map<unsigned long long,int>::iterator it;

    	// Print values
    	int count = 0;
    	for(auto const &it : ids) {
			iName = it.first;
    		i = it.second;

    		sprintf(buffer,"%d\t%llu\t%.12f",count,iName,M.get_object(i).e_x);
    		count++;
    		aux = string(buffer);
	    	write_line(filename,aux);
    	}
    }

    static void write_matrix(string filename, std::map<unsigned long long,int> &ids, Matrix2D<my_pvar_aux> &M) {
    	char buffer[100];
    	string aux;
    	int i;
    	unsigned long long iName;
    	std::map<unsigned long long,int>::iterator it;

    	// Print values
    	int count = 0;
    	for(const auto &it : ids) {
			iName = it.first;
    		i = it.second;

    		sprintf(buffer,"%d\t%llu\t",count,iName);
    		count++;
    		aux = string(buffer);
    		for(int k=0; k<M.get_size2()-1; k++) {
	    		sprintf(buffer,"%.12f\t",M.get_object(i,k).e_x);
    			aux += string(buffer);
    		}
    		sprintf(buffer,"%.12f",M.get_object(i,M.get_size2()-1).e_x);
    		aux += string(buffer);
	    	write_line(filename,aux);
    	}
    }

	static void write_matrix_users_sess(string filename, my_data &data, Matrix2D<my_pvar_aux> &M) {
    	char buffer[100];
    	int u;
    	int s;
    	int t;
    	string aux;

    	// Build auxiliary user_ids
    	std::vector<unsigned long long> user_ids = std::vector<unsigned long long>(data.Nusers);
    	for(auto const &it : data.user_ids) {
    		user_ids.at(it.second) = it.first;
    	}

    	// Build auxiliary session_ids
    	std::vector<unsigned long long> session_ids = std::vector<unsigned long long>(data.Nsessions);
    	for(auto const &it : data.session_ids) {
    		session_ids.at(it.second) = it.first;
    	}

    	// Print values
    	int count = 0;
    	for(const auto &it : data.trans_ids) {
    		u = it.first.u;
    		s = it.first.s;
    		t = it.second;

    		sprintf(buffer,"%d\t%llu\t%llu\t",count,user_ids.at(u),session_ids.at(s));
    		aux = string(buffer);
    		for(int k=0; k<M.get_size2()-1; k++) {
	    		sprintf(buffer,"%.12f\t",M.get_object(t,k).e_x);
    			aux += string(buffer);
    		}
    		sprintf(buffer,"%.12f",M.get_object(t,M.get_size2()-1).e_x);
    		aux += string(buffer);
	    	write_line(filename,aux);
	    	count++;
    	}
    }

    static void write_all_matrices(string label, my_data &data, const my_hyper &hyper, const my_param &param, my_pvar &pvar) {
	    // Print output matrices
        write_matrix(param.outdir+"/param_alpha"+label+".txt",data.item_ids,pvar.alpha);
	    if(param.flag_ppca) {
	        write_matrix_users_sess(param.outdir+"/param_rho"+label+".txt",data,pvar.rho);
	    } else {
	        write_matrix(param.outdir+"/param_rho"+label+".txt",data.item_ids,pvar.rho);
	        if(param.flag_userVec) {
		        write_matrix(param.outdir+"/param_theta"+label+".txt",data.user_ids,pvar.theta);
	        }
	        if(param.flag_itemIntercept) {
		        write_matrix(param.outdir+"/param_lambda0"+label+".txt",data.item_ids,pvar.lambda0);
	        }
	    }
	    if(param.flag_price>0) {
	    	write_matrix(param.outdir+"/param_gamma"+label+".txt",data.user_ids,pvar.gamma);
			write_matrix(param.outdir+"/param_beta"+label+".txt",data.item_ids,pvar.beta);
	    }
    }

    static void write_objective_function(const my_param &param, double logp) {
    	char buffer[100];
    	sprintf(buffer,"%d\t%.6f",param.it,logp);
    	string aux = string(buffer);
    	write_line(param.outdir+"/obj_function.txt",aux);
    }
};

#endif
