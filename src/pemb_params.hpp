#ifndef PEMB_PARAMS_HPP
#define PEMB_PARAMS_HPP

class my_param {
public:

	string datadir;
	string outdir;
	int K;
	unsigned long int seed;
	int rfreq;
	int saveCycle;
	int Niter;
	int negsamples;
	double zeroFactor;
	int batchsize;
	string label;
	int flag_avgContext;
	bool flag_itemIntercept;
	bool flag_ppca;
	bool flag_additiveMean;
	bool flag_gaussian;
	bool flag_bernoulli;
	bool flag_rmsprop;
	bool flag_adagrad;
	int flag_nsFreq;
	int flag_userVec;
	bool flag_binarizeContext;
	int flag_price;
	bool flag_day;
	int flag_tripEffects;
	bool flag_normPrice;
	bool flag_regularization;
	double eta;
	double gamma;
	double stdIni;
	double valTolerance;
	string iniPath;
	bool flag_iniThetaVal;
	double iniThetaVal;
	bool flag_iniPriceVal;
	double iniPriceVal;
	int valConsecutive;
	bool noVal;
	bool noTest;

	int lf_keepOnly;
	int lf_keepAbove;
	int lf_flag;

	int it;
	int n_val_decr;
	double prev_val_llh;
		
	my_param() {
		datadir = ".";
		outdir = ".";
		K = 50;
		seed = 0;
		rfreq = 10;
		saveCycle = 100;
		negsamples = 10;
		zeroFactor = 0.1;
		batchsize = -1;
		Niter = 1500;
		flag_avgContext = 1;
		flag_itemIntercept = false;
		flag_additiveMean = false;
		flag_ppca = false;
		flag_gaussian = false;
		flag_bernoulli = false;
		flag_rmsprop = false;
		flag_adagrad = true;
		flag_nsFreq = -1;
		flag_userVec = 0;
		flag_binarizeContext = false;
		flag_price = 0;
		flag_day = false;
		flag_normPrice = false;
		flag_regularization = true;
		flag_tripEffects = 0;
		label = "";
		eta = 0.1;
		gamma = 0.9;
		stdIni = -1.0;
		valTolerance = 0.000001;
		valConsecutive = 5;
		noVal = false;
		noTest = false;
		iniPath = "";
		flag_iniThetaVal = false;
		iniThetaVal = 0.0;
		flag_iniPriceVal = false;
		iniPriceVal = 0.0;

		lf_keepOnly = -1;
		lf_keepAbove = -1;
		lf_flag = 0;

		it = 0;
		prev_val_llh = 0.0;
		n_val_decr = 0;
	}
};

class my_hyper {
public:

	double s2rho;
	double s2alpha;
	double s2theta;
	double s2gamma;
	double s2beta;
	double s2trip;
	double mean_gamma;
	double mean_beta;
	double mean_trip;
	double s2noise;
	double mean_xday;
	double s2xday;

	my_hyper() {
		s2rho = 1.0;
		s2alpha = 1.0;
		s2theta = 1.0;
		s2noise = 1.0;
		s2gamma = 1.0;
		s2beta = 1.0;
		s2trip = 1.0;
		mean_gamma = -5.0;
		mean_beta = -5.0;
		s2xday = 0.1;
		mean_xday = -5.0;
		mean_trip = -5.0;
	}
};

class my_data_aux {
public:
	unsigned int T;			// Number of (user,session,item) triplets
	unsigned int *y_user;	// User idx per (user,session,item) triplet
	unsigned int *y_item;	// Item idx per (user,session,item) triplet
	unsigned int *y_sess;	// Session idx per (user,session,item) triplet
	unsigned int *y_rating;	// Value (rating or #units) per (user,session,item) triplet
	unsigned int *y_trans;	// Transaction idx per (user,session,item) triplet

	my_data_aux() {
		T = 0;
		y_user = nullptr;
		y_item = nullptr;
		y_sess = nullptr;
		y_rating = nullptr;
		y_trans = nullptr;
	}

	my_data_aux(unsigned int T_) {
		allocate_all(T_);
	}

	~my_data_aux() {
		delete_all();
	}

	inline void allocate_all(unsigned int T_) {
		T = T_;
		y_user = new unsigned int[T];
		y_item = new unsigned int[T];
		y_sess = new unsigned int[T];
		y_rating = new unsigned int[T];
		y_trans = new unsigned int[T];
	}

	inline void delete_all() {
		delete [] y_user;
		delete [] y_item;
		delete [] y_sess;
		delete [] y_rating;
		delete [] y_trans;
	}

	my_data_aux & operator=(const my_data_aux &rhs) {
    	// Check for self-assignment!
	    if(this!=&rhs) {
	    	// deallocate memory
			delete_all();
			// allocate memory for the contents of rhs
			allocate_all(rhs.T);
			// copy values from rhs
			for(unsigned int n=0; n<rhs.T; n++) {
				y_user[n] = rhs.y_user[n];
				y_item[n] = rhs.y_item[n];
				y_sess[n] = rhs.y_sess[n];
				y_rating[n] = rhs.y_rating[n];
				y_trans[n] = rhs.y_trans[n];
			}
		}
		return *this;
	}
};

class hpf_trans_aux {
public:
	int u;
	int s;

	hpf_trans_aux & operator=(const hpf_trans_aux &rhs) {
    	// Check for self-assignment!
	    if (this!=&rhs) {
	    	u = rhs.u;
	    	s = rhs.s;
		}
		return *this;
	}

	hpf_trans_aux(int uu, int ss) {
		u = uu;
		s = ss;
	}

	bool operator==(const hpf_trans_aux &b) const {
		return (u==b.u)&&(s==b.s);
	}

	bool operator!=(const hpf_trans_aux &b) const {
		return (u!=b.u)||(s!=b.s);
	}

	bool operator<(const hpf_trans_aux &b) const {
		return (u<b.u)||(u==b.u && s<b.s);
	}

	bool operator<=(const hpf_trans_aux &b) const {
		return (u<b.u)||(u==b.u && s<=b.s);
	}

	bool operator>(const hpf_trans_aux &b) const {
		return (u>b.u)||(u==b.u && s>b.s);
	}

	bool operator>=(const hpf_trans_aux &b) const {
		return (u>b.u)||(u==b.u && s>=b.s);
	}
};

class my_data {
public:

	// Integers
	int Nitems;		// Number of items
	int Nusers;		// Number of users
	int Nsessions;	// Number of sessions
	int Ntrans;		// Number of transactions (a transaction is a (user,session) pair)
	int Ndays;		// Number of calendar days
	int Nweekdays;  // Number of weekdays (typically 7)
	int NuserGroups;	// Number of user groups
	int NitemGroups;	// Number of item groups

	// Structs with the actual data
	my_data_aux obs;			// Observations (train)
	my_data_aux obs_test;		// Observations (test)
	my_data_aux obs_val;		// Observations (validation)

	// Prices for each session and item
	Matrix2D<double> price_is;	  // Prices
	Matrix1D<double> price_avg_i; // Average price for each item (across all sessions)

	// sum(log(y!))
	double sum_log_yfact;
	Matrix1D<double> sum_log_yfact_per_trans;

	// Mapping from ids to indices
	std::map<unsigned long long, int> item_ids;		// Map containing the item id's
	std::map<unsigned long long, int> user_ids;		// Map containing the user id's
	std::map<unsigned long long, int> session_ids;	// Map containing the session id's
	std::map<hpf_trans_aux,int> trans_ids;			// Map containing the transaction id's (for train+validation)
	std::map<unsigned long long,int> day_ids;			// Map containing the day id's
	std::map<unsigned long long,int> weekday_ids;		// Map containing the weekday id's
	std::map<unsigned long long,int> usergroup_ids;		// Map containing the usergroups id's
	std::map<unsigned long long,int> itemgroup_ids;		// Map containing the itemgroups id's

	// Misc lists
	Matrix1D<std::vector<int>> sessions_per_user;	// For each user, list of sessions in which that user appears
	Matrix1D<std::vector<int>> trans_per_item;		// For each item, list of transactions in which that item appears
	Matrix1D<std::vector<int>> trans_per_user;		// For each user, list of transactions in which that user appears
	Matrix1D<int> user_per_trans;					// For each transaction, user to which it corresponds
	Matrix1D<int> session_per_trans;				// For each transaction, session to which it corresponds
	Matrix1D<std::vector<int>> lines_per_trans;		// For each transaction, list of "lines" in train.tsv referring to that transaction
	Matrix1D<std::vector<int>> items_per_trans;		// For each transaction, list of items in that transaction
	int maxNi;										// Maximum #items of all transactions
	Matrix1D<std::vector<int>> lines_per_item;		// For each item, list of "lines" in train.tsv in which that item appears
	Matrix1D<std::vector<int>> lines_per_user;		// For each user, list of "lines" in train.tsv in which that user appears
	Matrix1D<int> Nitems_per_trans;					// For each transaction, sum(y)
	Matrix1D<int> sum_sizetrans_per_item;			// For each item, sum_t(length(t)) for all the transactions t in which that item appears
	Matrix1D<int> group_per_user;					// For each user, group to which she belongs
	Matrix1D<int> group_per_item;					// For each item, group to which it belongs
	Matrix1D<std::vector<int>> users_per_group;		// For each usergroup, list of the users that belong to that group
	Matrix1D<std::vector<int>> items_per_group;		// For each itemgroup, list of the items that belong to that group
	Matrix1D<int> day_per_session;					// For each session, day_id to which it corresponds
	Matrix1D<int> weekday_per_session;				// For each session, weekday_id to which it corresponds
	Matrix1D<double> hour_per_session;				// For each session, hour of the day
	Matrix1D<std::vector<int>> sessions_per_day;	// For each day, a list of sessions
	Matrix1D<std::vector<int>> sessions_per_weekday;// For each weekday, a list of sessions
	Matrix3D<std::vector<int>> lines_per_xday;		// For each day, a list of "lines" in train.tsv

	// Misc auxiliary counts
	Matrix2D<int> count_trans_per_xday;				// Each entry = [#trips that all users from group g_u make in day d]

	// Valid list of items
	std::vector<unsigned long long> valid_items;	// Used to remove low-frequency items

	// Lists for transactions in test set (same definition as above, but for test set only)
	int test_Ntrans;
	Matrix1D<std::vector<int>> test_lines_per_trans;
	Matrix1D<std::vector<int>> test_items_per_trans;
	Matrix1D<int> test_Nitems_per_trans;
	Matrix1D<int> test_user_per_trans;
	Matrix1D<int> test_session_per_trans;
	Matrix1D<std::vector<int>> test_sessions_per_user;
	std::map<hpf_trans_aux,int> test_trans_ids;			// Map containing the transaction id's (for test)

	// Frequency and counts of items
	double *uniform_dist;
	double *unigram_dist;
	double *unigram_dist_power;
	gsl_ran_discrete_t *negsampling_dis;

	my_data() {
		Nitems = 0;
		Nusers = 0;
		Nsessions = 0;
		Ntrans = 0;
		Ndays = 0;
		Nweekdays = 0;
		NuserGroups = 0;
		NitemGroups = 0;
		test_Ntrans = 0;
		maxNi = 0;
		uniform_dist = nullptr;
		unigram_dist = nullptr;
		unigram_dist_power = nullptr;
		negsampling_dis = nullptr;
	}

	~my_data() {
		delete [] uniform_dist;
		delete [] unigram_dist;
		delete [] unigram_dist_power;
		gsl_ran_discrete_free(negsampling_dis);
	}

	inline double get_price(int i, int s, const my_param &param) {
		return((param.flag_price>0)?price_is.get_object(i,s):0.0);
	}

	inline int get_order_in_trans(int i, int t) {
		bool found = false;
		int y = -1;
		int j;
		int count = 0;
		std::vector<int>::iterator iter = items_per_trans.get_object(t).begin();
		while(!found && iter!=items_per_trans.get_object(t).end()) {
			j = *iter;
			if(j==i) {
				found = true;
				y = count;
			}
			count++;
			++iter;
		}
		return y;
	} 

	inline double get_items_in_trans(int i, int t) {
		bool found = false;
		int ll;
		double y = -1.0;
		std::vector<int>::iterator iter = lines_per_trans.get_object(t).begin();
		while(!found && iter!=lines_per_trans.get_object(t).end()) {
			ll = *iter;
			if(i==static_cast<int>(obs.y_item[ll])) {
				found = true;
				y = static_cast<double>(obs.y_rating[ll]);
			}
			++iter;
		}
		return y;
	}

	void compute_unigram_distributions(const my_param &param) {
		uniform_dist = new double[Nitems];
		unigram_dist = new double[Nitems];
		unigram_dist_power = new double[Nitems];

		// Count item occurrences
		for(unsigned int t=0; t<obs.T; t++) {
			int i = obs.y_item[t];
			unigram_dist[i] += 1.0;
		}
		// Normalize counts
		for(int i=0; i<Nitems; i++) {
			unigram_dist[i] /= static_cast<double>(obs.T);
			uniform_dist[i] = 1.0/static_cast<double>(Nitems);
		}
		// Raise to the power of 3/4
		double suma = 0.0;
		for(int i=0; i<Nitems; i++) {
			unigram_dist_power[i] = pow(unigram_dist[i],0.75);
			suma += unigram_dist_power[i];
		}
		// Normalize the power'ed counts
		for(int i=0; i<Nitems; i++) {
			unigram_dist_power[i] /= suma;
		}
		// Create preprocessed gsl_ran_discrete_t variables for faster negative sampling
		if(param.flag_nsFreq==-1) {
			negsampling_dis = gsl_ran_discrete_preproc(Nitems,uniform_dist);
		} else if(param.flag_nsFreq==0) {
			negsampling_dis = gsl_ran_discrete_preproc(Nitems,unigram_dist);
		} else if(param.flag_nsFreq==1) {
			negsampling_dis = gsl_ran_discrete_preproc(Nitems,unigram_dist_power);
		} else {
			std::cerr << "[ERR] Wrong value for -nsFreq" << endl;
			assert(0);
		}
	}

	void create_transactions_train(const my_param &param) {
		int u;
		int s;
		Ntrans = 0;
		for(unsigned int t=0; t<obs.T; t++) {
			u = obs.y_user[t];
			s = obs.y_sess[t];
			hpf_trans_aux tid = hpf_trans_aux(u,s);
			if(trans_ids.find(tid)==trans_ids.end()) {
				trans_ids.insert(std::pair<hpf_trans_aux,int>(tid,Ntrans));
				Ntrans += 1;
			}
		}
	}

	void create_sessions_per_user(const my_param &param) {
		sessions_per_user = Matrix1D<std::vector<int>>(Nusers);
		int u;
		int s;
		for(auto const &it : trans_ids) {
			u = it.first.u;
			s = it.first.s;
			sessions_per_user.get_object(u).push_back(s);
		}
	}

	void create_other_data_structs(const my_param &param) {
		int u;
		int i;
		int idx_trans;
		int s;

		// Compute sum_log_yfact
		sum_log_yfact = 0.0;
		for(unsigned int t=0; t<obs.T; t++) {
			sum_log_yfact += my_logfactorial(obs.y_rating[t]);
		}

		// Create session_per_trans, user_per_trans
		session_per_trans = Matrix1D<int>(Ntrans);
		user_per_trans = Matrix1D<int>(Ntrans);
		for(auto const &it : trans_ids) {
			u = it.first.u;
			s = it.first.s;
			idx_trans = it.second;
			user_per_trans.set_object(idx_trans,u);
			session_per_trans.set_object(idx_trans,s);
		}

		// Create lists mapping from lines<-->trans, and create mapping item-->transactions
		trans_per_item = Matrix1D<std::vector<int>>(Nitems);
		lines_per_item = Matrix1D<std::vector<int>>(Nitems);
		lines_per_user = Matrix1D<std::vector<int>>(Nusers);
		lines_per_trans = Matrix1D<std::vector<int>>(Ntrans);
		items_per_trans = Matrix1D<std::vector<int>>(Ntrans);
		Nitems_per_trans = Matrix1D<int>(Ntrans);
		set_to_zero(Nitems_per_trans);
		for(unsigned int t=0; t<obs.T; t++) {
			u = obs.y_user[t];
			s = obs.y_sess[t];
			i = obs.y_item[t];

			idx_trans = get_transaction(u,s);
			trans_per_item.get_object(i).push_back(idx_trans);
			lines_per_trans.get_object(idx_trans).push_back(t);
			items_per_trans.get_object(idx_trans).push_back(i);
			lines_per_item.get_object(i).push_back(t);
			lines_per_user.get_object(u).push_back(t);
			obs.y_trans[t] = idx_trans;

			int aux_int = Nitems_per_trans.get_object(idx_trans);
			aux_int += obs.y_rating[t];
			Nitems_per_trans.set_object(idx_trans,aux_int);
		}

		// Compute sum_log_yfact_per_trans
		sum_log_yfact_per_trans = Matrix1D<double>(Ntrans);
		for(int t=0; t<Ntrans; t++) {
			double sum = 0.0;
			for(int &ll : lines_per_trans.get_object(t)) {
				sum += my_logfactorial(obs.y_rating[ll]);
			}
			sum_log_yfact_per_trans.set_object(t,sum);
		}

		// Set maxNi
		for(int t=0; t<Ntrans; t++) {
			int Nt = items_per_trans.get_object(t).size();
			maxNi = my_max(Nt,maxNi);
		}

		// Create obs_val.y_trans[]
		if(!param.noVal){
			for(unsigned int t=0; t<obs_val.T; t++) {
				u = obs_val.y_user[t];
				s = obs_val.y_sess[t];
				i = obs_val.y_item[t];

				idx_trans = get_transaction(u,s);
				if(idx_trans<0) {
					std::cerr << "[WARN] Line " << (t+1) << " of validation.tsv contains a (user,session) pair that is not present in train.tsv." \
					          << " This line will be ignored" << endl;
				}
				obs_val.y_trans[t] = idx_trans;
			}
		}

		// Create trans_per_user
		trans_per_user = Matrix1D<std::vector<int>>(Nusers);
		for(int t=0; t<Ntrans; t++) {
			u = user_per_trans.get_object(t);
			trans_per_user.get_object(u).push_back(t);
		}

		// Create sum_sizetrans_per_item
		sum_sizetrans_per_item = Matrix1D<int>(Nitems);
		for(int i=0; i<Nitems; i++) {
			int cc = 0;
			// for all transactions in which item i appears
			for(int t : trans_per_item.get_object(i)) {
				cc += items_per_trans.get_object(t).size();
			}
			sum_sizetrans_per_item.set_object(i,cc);
		}
	}

	void create_count_trans_per_xday(const my_param &param) {
		// Create count_trans_per_xday
		if(param.flag_day) {
			count_trans_per_xday = Matrix2D<int>(NuserGroups,Ndays);
			set_to_zero(count_trans_per_xday);
			// iterate over all transactions
			for(auto const &iter : trans_ids) {
				int u = iter.first.u;
				int s = iter.first.s;
				int d = day_per_session.get_object(s);
				int g_u = group_per_user.get_object(u);
				int val = 1+count_trans_per_xday.get_object(g_u,d);
				count_trans_per_xday.set_object(g_u,d,val);
			}
		}
	}

	void create_transactions_test(const my_param &param) {
		if(param.noTest) {
			return;
		}

		int u;
		int s;
		int i;
		int aux;

		// Create test_trans_ids
		test_Ntrans = 0;
		for(unsigned int t=0; t<obs_test.T; t++) {
			u = obs_test.y_user[t];
			s = obs_test.y_sess[t];
			hpf_trans_aux tid = hpf_trans_aux(u,s);
			if(test_trans_ids.find(tid)==test_trans_ids.end()) {
				test_trans_ids.insert(std::pair<hpf_trans_aux,int>(tid,test_Ntrans));
				test_Ntrans += 1;
			}
		}

		// Create test_sessions_per_user
		test_sessions_per_user = Matrix1D<std::vector<int>>(Nusers);
		for(auto const &it : test_trans_ids) {
			u = it.first.u;
			s = it.first.s;
			test_sessions_per_user.get_object(u).push_back(s);
		}

		// Create test_user_per_trans and test_session_per_trans
		test_user_per_trans = Matrix1D<int>(test_Ntrans);
		test_session_per_trans = Matrix1D<int>(test_Ntrans);
		for(auto const &it : test_trans_ids) {
			u = it.first.u;
			s = it.first.s;
			int t = it.second;
			test_user_per_trans.set_object(t,u);
			test_session_per_trans.set_object(t,s);
		}

		// Create test_lines_per_trans, test_items_per_trans, and obs_test.y_trans
		test_lines_per_trans = Matrix1D<std::vector<int>>(test_Ntrans);
		test_items_per_trans = Matrix1D<std::vector<int>>(test_Ntrans);
		test_Nitems_per_trans = Matrix1D<int>(test_Ntrans);
		set_to_zero(test_Nitems_per_trans);
		for(unsigned int ll=0; ll<obs_test.T; ll++) {
			u = obs_test.y_user[ll];
			s = obs_test.y_sess[ll];
			i = obs_test.y_item[ll];
			int t = test_get_transaction(u,s);
			if(t<0) {
				std::cerr << "[ERR] Transaction not found (my_data::create_transactions_test)" << endl;
				assert(0);
			}
			obs_test.y_trans[ll] = t;
			test_lines_per_trans.get_object(t).push_back(ll);
			test_items_per_trans.get_object(t).push_back(i);
			aux = test_Nitems_per_trans.get_object(t);
			aux += static_cast<int>(obs_test.y_rating[ll]);
			test_Nitems_per_trans.set_object(t,aux);
		}
	}

	inline int get_transaction(int u, int s) {
		int t = -1;
		hpf_trans_aux tid = hpf_trans_aux(u,s);
		std::map<hpf_trans_aux,int>::const_iterator iter = trans_ids.find(tid);
		if(iter!=trans_ids.end()) {
			t = iter->second;
		}
		return t;
	}

	inline int test_get_transaction(int u, int s) {
		int t = -1;
		hpf_trans_aux tid = hpf_trans_aux(u,s);
		std::map<hpf_trans_aux,int>::const_iterator iter = test_trans_ids.find(tid);
		if(iter!=test_trans_ids.end()) {
			t = iter->second;
		}
		return t;
	}

	inline void set_to_zero(Matrix1D<int> &m) {
		for(int n=0; n<m.get_size1(); n++) {
			m.set_object(n,0);
		}
	}

	inline void set_to_zero(Matrix2D<int> &m) {
		for(int n=0; n<m.get_size1(); n++) {
			for(int k=0; k<m.get_size2(); k++) {
				m.set_object(n,k,0);
			}
		}
	}
};

class my_pvar_aux {
public:
	double e_x;
	double e_log;
	double grad;
	double grad_neg;
	double Gt;
	bool flag_positive;

	my_pvar_aux() {
		e_x = 0.0;
		e_log = 0.0;
		grad = 0.0;
		grad_neg = 0.0;
		Gt = 0.0;
		flag_positive = false;
	}

	void initialize_value(double vv, bool pp) {
		flag_positive = pp;
		if(flag_positive) {
			e_x = vv;
	        e_log = my_log(vv);
		} else {
			e_x = vv;
		}
	}

	void initialize_random(gsl_rng *semilla, double vv, double stdIni, bool pp) {
		flag_positive = pp;
		double aux;
		if(flag_positive) {
			double shape = my_pow2(vv)/my_pow2(stdIni);
			double rate = shape/vv;			
			aux = gsl_ran_gamma(semilla,shape,1.0/rate);
	        aux = (aux<1e-5?1e-5:aux);
	        e_x = aux;
	        e_log = my_log(aux);
		} else {
			e_x = vv+gsl_ran_gaussian(semilla,stdIni);
		}
	}

	inline void set_grad_to_zero() {
		grad = 0.0;
		grad_neg = 0.0;
	}

	inline double set_grad_to_prior(double ss2) {
		double logp = set_grad_to_prior(0.0,ss2);
		return logp;
	}

	inline double set_grad_to_prior(double mm, double ss2) {
		double logp;
		if(flag_positive) {
			grad = -(e_log-mm)/ss2;
			logp = 0.5*(e_log-mm)*grad;
		} else {
			grad = -(e_x-mm)/ss2;
			logp = 0.5*(e_x-mm)*grad;
		}
		grad_neg = 0.0;
		return logp;
	}

	inline void increase_grad(double val) {
		grad += val;
	}

	inline void increase_grad_neg(double val) {
		grad_neg += val;
	}

	inline void scale_add_grad_neg(double sviFactor) {
		grad_neg *= sviFactor;
		grad += grad_neg;
	}

	inline void take_grad_step(const my_param &param, double stepsize, double gamma, int it) {
		// Convert grad into gradient of the log
		if(flag_positive) {
			grad *= e_x;
		}
		// Increase Gt
		if(param.flag_adagrad) {
			Gt += my_pow2(grad);
		} else if(param.flag_rmsprop) {
			gamma = (it==0)?0.0:gamma;
			Gt = (1.0-gamma)*my_pow2(grad)+gamma*Gt;
		} else {
			std::cerr << "[ERR] Choose either rmsprop or adagrad" << endl;
			assert(0);
		}
		// Take gradient step
		if(flag_positive) {
			if(Gt>0.0) {
				e_log += stepsize*grad/sqrt(Gt);
				e_log = (e_log<-18.4207?-18.4207:e_log);
				e_log = (e_log>3.9120?3.9120:e_log);
				e_x = my_exp(e_log);
			}
		} else {
			if(Gt>0.0) {
				e_x += stepsize*grad/sqrt(Gt);
				// e_log is not updated because it isn't used
			}
		}
	}
};

class my_pvar {
public:
	Matrix2D<my_pvar_aux> rho;      // embedding vectors
	Matrix2D<my_pvar_aux> alpha;    // context vectors
	Matrix2D<my_pvar_aux> theta;    // user vectors
	Matrix1D<my_pvar_aux> lambda0;  // item intercepts
	Matrix2D<double> sum_alpha;     // auxiliary variable (sum of alpha's in the context)
	Matrix2D<my_pvar_aux> gamma;    // price sensitivity vectors (per-user)
	Matrix2D<my_pvar_aux> beta;		// price sensitivity vectors (per-item)
	Matrix3D<my_pvar_aux> x_gid;	// per-day effects
	Matrix2D<my_pvar_aux> theta_trip;	// per-trip variables
	Matrix2D<my_pvar_aux> beta_trip;	// per-item variables corresponding to per-trip effects

	my_pvar(const my_data &data, const my_param &param) {
		int sizeRho = (param.flag_ppca?data.Ntrans:data.Nitems);
		int sizeSum = (param.flag_ppca?data.maxNi:param.K);
		
		rho = Matrix2D<my_pvar_aux>(sizeRho,param.K);
		alpha = Matrix2D<my_pvar_aux>(data.Nitems,param.K);
		sum_alpha = Matrix2D<double>(data.Ntrans,sizeSum);
		if(param.flag_userVec>0) {
			theta = Matrix2D<my_pvar_aux>(data.Nusers,param.K);
		}
		if(param.flag_itemIntercept) {
			lambda0 = Matrix1D<my_pvar_aux>(data.Nitems);
		}
		if(param.flag_price>0) {
			gamma = Matrix2D<my_pvar_aux>(data.Nusers,param.flag_price);
			beta = Matrix2D<my_pvar_aux>(data.Nitems,param.flag_price);
		}
		if(param.flag_day) {
			x_gid = Matrix3D<my_pvar_aux>(data.NuserGroups,data.NitemGroups,data.Ndays);
		}
		if(param.flag_tripEffects>0) {
			theta_trip = Matrix2D<my_pvar_aux>(data.Ntrans,param.flag_tripEffects);
			beta_trip = Matrix2D<my_pvar_aux>(data.Nitems,param.flag_tripEffects);
		}
	}

	void initialize_all(gsl_rng *semilla, my_data &data, const my_param &param, const my_hyper &hyper) {
		double val_ini = param.flag_additiveMean?1.0/sqrt(param.K):0.0;
		// initialize rho
		if(param.iniPath=="") {
			// initialize randomly
			initialize_matrix_randomly(semilla,rho,val_ini,param.stdIni,param.flag_additiveMean);
		} else {
			// initialize from file
			if(!param.flag_ppca) {
			  	initialize_from_file(param.iniPath+"/param_rho.txt",data.item_ids,rho,param,param.flag_additiveMean);
			} else {
			  	initialize_from_file_ppca(param.iniPath+"/param_rho.txt",data,rho,param,false);
			}
		}

		// initialize alpha
		if(param.iniPath=="") {
			// initialize randomly
			initialize_matrix_randomly(semilla,alpha,val_ini,param.stdIni,param.flag_additiveMean);
		} else {
			// initialize from file
			initialize_from_file(param.iniPath+"/param_alpha.txt",data.item_ids,alpha,param,param.flag_additiveMean);
		}

		// initialize theta
		if(param.flag_userVec>0) {
			if(param.flag_iniThetaVal) {
				// initialize to fixed value
				for(int u=0; u<data.Nusers; u++) {
					for(int k=0; k<param.K; k++) {
						theta.get_object(u,k).initialize_value(param.iniThetaVal,param.flag_additiveMean);
					}
				}
			} else if(param.iniPath=="") {
				// initialize randomly
				initialize_matrix_randomly(semilla,theta,val_ini,param.stdIni,param.flag_additiveMean);
			} else {
				// initialize from file
				initialize_from_file(param.iniPath+"/param_theta.txt",data.user_ids,theta,param,param.flag_additiveMean);
			}
		}
		
		// initialize intercepts lambda0
		if(param.flag_itemIntercept) {
			if(param.iniPath=="") {
				// initialize randomly
				initialize_matrix_randomly(semilla,lambda0,val_ini,param.stdIni,param.flag_additiveMean);
			} else {
				// initialize from file
				initialize_from_file(param.iniPath+"/param_lambda0.txt",data.item_ids,lambda0,param,param.flag_additiveMean);
			}
		}

		// initialize price sensitivities (beta, gamma)
		if(param.flag_price>0) {
			if(param.flag_iniPriceVal) {
				// initialize to fixed value
				for(int u=0; u<data.Nusers; u++) {
					for(int k=0; k<param.flag_price; k++) {
						gamma.get_object(u,k).initialize_value(param.iniPriceVal,true);
					}
				}
				for(int i=0; i<data.Nitems; i++) {
					for(int k=0; k<param.flag_price; k++) {
						beta.get_object(i,k).initialize_value(param.iniPriceVal,true);
					}
				}
			} else if(param.iniPath=="") {
				// initialize randomly
				initialize_matrix_randomly(semilla,gamma,0.1/sqrt(param.flag_price),param.stdIni,true);
				initialize_matrix_randomly(semilla,beta,0.1/sqrt(param.flag_price),param.stdIni,true);
			} else {
				// initialize from file
				initialize_from_file(param.iniPath+"/param_gamma.txt",data.user_ids,gamma,param,true);
				initialize_from_file(param.iniPath+"/param_beta.txt",data.item_ids,beta,param,true);
			}

			// initialize per-day effect variables
			if(param.flag_day) {
				if(param.iniPath=="") {
					// initialize randomly
					initialize_matrix_randomly(semilla,x_gid,val_ini,param.stdIni,param.flag_additiveMean);
				} else {
					// initialize from file
					initialize_from_file_xday(param.iniPath+"/param_xday.txt",data.usergroup_ids,data.itemgroup_ids,data.day_ids,x_gid,param,param.flag_additiveMean);
				}
			}
		}

		// initialize trip effects
		if(param.flag_tripEffects>0) {
			val_ini = 1.0/sqrt(static_cast<double>(param.flag_tripEffects));
			if(param.iniPath=="") {
				// initialize randomly
				initialize_matrix_randomly(semilla,theta_trip,val_ini,param.stdIni,true);
				initialize_matrix_randomly(semilla,beta_trip,val_ini,param.stdIni,true);
			} else {
				std::cerr << "[WARN] Trip effects will be initialized randomly" << endl;
			}
		}
	}

	static void initialize_matrix_randomly(gsl_rng *semilla, Matrix1D<my_pvar_aux> &M, double vv, double ss, bool pp) {
		for(int i=0; i<M.get_size1(); i++) {
			M.get_object(i).initialize_random(semilla,vv,ss,pp);
		}
	}

	static void initialize_matrix_randomly(gsl_rng *semilla, Matrix2D<my_pvar_aux> &M, double vv, double ss, bool pp) {
		for(int i=0; i<M.get_size1(); i++) {
			for(int k=0; k<M.get_size2(); k++) {
				M.get_object(i,k).initialize_random(semilla,vv,ss,pp);
			}
		}
	}

	static void initialize_matrix_randomly(gsl_rng *semilla, Matrix3D<my_pvar_aux> &M, double vv, double ss, bool pp) {
		for(int i=0; i<M.get_size1(); i++) {
			for(int k=0; k<M.get_size2(); k++) {
				for(int m=0; m<M.get_size3(); m++) {
					M.get_object(i,k,m).initialize_random(semilla,vv,ss,pp);
				}
			}
		}
	}

	static void initialize_from_file(string fname, const std::map<unsigned long long, int> &ids, Matrix1D<my_pvar_aux> &M, const my_param &param, bool pp) {
		FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
	  	unsigned long long id;
	  	int ll;
	  	int idx;
	  	double vv;
		while(!feof(fin)){
			// Read a line
			fscanf(fin,"%d\t%llu\t",&ll,&id);
			std::map<unsigned long long,int>::const_iterator iter = ids.find(id);
			if(iter==ids.end()) {
		  		std::cerr << "[ERR] Error reading line " << ll << " of " << fname << endl;
		  		std::cerr << "      Index not found in data files" << endl;
		  		assert(0);
			}
			idx = iter->second;
			// Read the value
			fscanf(fin,"%lf\n",&vv);
			M.get_object(idx).initialize_value(vv,pp);
		}
	  	fclose(fin);
	}

	static void initialize_from_file(string fname, const std::map<unsigned long long, int> &ids, Matrix2D<my_pvar_aux> &M, const my_param &param, bool pp) {
		FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
	  	unsigned long long id;
	  	int ll;
	  	int idx;
	  	double vv;
		while(!feof(fin)){
			// Read a line
			fscanf(fin,"%d\t%llu\t",&ll,&id);
			std::map<unsigned long long,int>::const_iterator iter = ids.find(id);
			if(iter==ids.end()) {
		  		std::cerr << "[ERR] Error reading line " << ll << " of " << fname << endl;
		  		std::cerr << "      Index not found in data files" << endl;
		  		assert(0);
			}
			idx = iter->second;
			for(int k=0; k<M.get_size2(); k++) {
				// Read a value
				if(k<M.get_size2()-1) {
					fscanf(fin,"%lf\t",&vv);
				} else {
					fscanf(fin,"%lf\n",&vv);
				}
				M.get_object(idx,k).initialize_value(vv,pp);
			}

		}
	  	fclose(fin);
	}

	static void initialize_from_file_xday(string fname, const std::map<unsigned long long, int> &usergroup_ids, const std::map<unsigned long long, int> &itemgroup_ids, const std::map<unsigned long long, int> &day_ids, \
										  Matrix3D<my_pvar_aux> &M, const my_param &param, bool pp) {
		FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
	  	unsigned long long gu_id;
	  	unsigned long long gi_id;
	  	unsigned long long d_id;
	  	int ll;
	  	int g_u;
	  	int g_i;
	  	int d;
	  	double vv;
		while(!feof(fin)){
			// Read a line
			fscanf(fin,"%d\t%llu\t%llu\t%llu\t%lf\n",&ll,&gu_id,&gi_id,&d_id,&vv);
			std::map<unsigned long long,int>::const_iterator iter_gu = usergroup_ids.find(gu_id);
			std::map<unsigned long long,int>::const_iterator iter_gi = itemgroup_ids.find(gi_id);
			std::map<unsigned long long,int>::const_iterator iter_d = day_ids.find(d_id);
			if(iter_gu==usergroup_ids.end() || iter_gi==itemgroup_ids.end() || iter_d==day_ids.end()) {
		  		std::cerr << "[ERR] Error reading line " << ll << " of " << fname << endl;
		  		std::cerr << "      One of the indices was not found in data files" << endl;
		  		assert(0);
			}
			g_u = iter_gu->second;
			g_i = iter_gi->second;
			d = iter_d->second;
			M.get_object(g_u,g_i,d).initialize_value(vv,pp);
		}
	  	fclose(fin);
	}	

	static void initialize_from_file_ppca(string fname, my_data &data, Matrix2D<my_pvar_aux> &M, const my_param &param, bool pp) {
		FILE *fin = fopen(fname.c_str(),"r");
	  	if(!fin) {
	  		std::cerr << "[ERR] Unable to open " << fname << endl;
	  		assert(0);
	  	}
	  	unsigned long long id_u;
	  	unsigned long long id_s;
	  	int ll;
	  	int idx_u;
	  	int idx_s;
	  	int t;
	  	double vv;
		while(!feof(fin)){
			// Read (user,session) pair
			fscanf(fin,"%d\t%llu\t%llu\t",&ll,&id_u,&id_s);
			std::map<unsigned long long,int>::const_iterator iter_u = data.user_ids.find(id_u);
			std::map<unsigned long long,int>::const_iterator iter_s = data.session_ids.find(id_s);
			if(iter_u==data.user_ids.end() || iter_s==data.session_ids.end()) {
		  		std::cerr << "[ERR] Error reading line " << ll << " of " << fname << endl;
		  		std::cerr << "      Index not found in data files" << endl;
		  		assert(0);
			}
			idx_u = iter_u->second;
			idx_s = iter_s->second;
			// Map (user,session) pair to transaction index t
			t = data.get_transaction(idx_u,idx_s);
			// Sanity check (transaction must be found)
			if(t<0) {
		  		std::cerr << "[ERR] Error reading line " << ll << " of " << fname << endl;
		  		std::cerr << "      (user,session) pair not found in data files" << endl;
		  		assert(0);
			}
			// Once mapped, read values
			for(int k=0; k<M.get_size2(); k++) {
				// Read a value
				if(k<M.get_size2()-1) {
					fscanf(fin,"%lf\t",&vv);
				} else {
					fscanf(fin,"%lf\n",&vv);
				}
				M.get_object(t,k).initialize_value(vv,pp);
			}
		}
	  	fclose(fin);
	}
};

#endif
