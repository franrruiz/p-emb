#include "my_headers_pemb.hpp"
#include "my_gsl_utilities.hpp"
#include "my_matrices.hpp"
#include "pemb_params.hpp"
#include "pemb_io.hpp"
#include "pemb_inference.hpp"
//#include "pemb_inference_ppca.hpp"

int main (int argc, char *argv[]) {

    my_param param;
    my_hyper hyper;
    my_data data;

    /******** Read options from input (argc, argv) ********/
    std::cout << "Initializing program..." << endl;
    my_input::read_input_from_command_line(argc,argv,data,param,hyper);

    /******** Read data from file ********/
    std::cout << "Reading data..." << endl;
    // Read train.tsv, validation.tsv and train.tsv
    my_input::read_data_file(data,param);
    // Create auxiliary data structures
    std::cout << "Computing auxiliary data structures..." << endl;
    data.create_transactions_train(param);
    data.create_sessions_per_user(param);
    data.create_other_data_structs(param);
    data.compute_unigram_distributions(param);
    std::cout << "Computing auxiliary data structures for test set..." << endl;
    data.create_transactions_test(param);
    // Create output folder
    my_output::create_output_folder(data,param);
    // Print vocabulary
    my_output::write_vocab(data,param);
    // Read prices
    if(param.flag_price>0) {
        std::cout << "Reading prices..." << endl;
        my_input::read_price_file(data,param);
        my_input::normalize_prices(data,param);
        my_output::write_normalized_prices(data,param);
        my_input::take_log_prices(data,param);
    }
    // Read usergroups,itemgroups,sess_days
    if(param.flag_day) {
        // read userGroup.tsv
        std::cout << "Reading user groups..." << endl;
        my_input::read_usergroups(data,param);
        // read itemGroup.tsv
        std::cout << "Reading item groups..." << endl;
        my_input::read_itemgroups(data,param);
        // read sess_days.tsv
        std::cout << "Reading session/days mapping..." << endl;
        my_input::read_sess_days(data,param);
        // Create lines_per_xday and count_trans_per_xday
        my_input::create_lines_per_xday(data,param);
        data.create_count_trans_per_xday(param);
    }

    /******** Write log file ********/
    my_output::create_log_file(data,param,hyper);

    /******** Set the seeds ********/
    srand48(param.seed);
    gsl_rng *semilla = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(semilla,param.seed);

    /******** Allocate memory for variational parameters ********/
    std::cout << "Initializing latent parameters..." << endl;
    my_pvar pvar(data,param);
    // Initialize variational parameters
    pvar.initialize_all(semilla,data,param,hyper);

    /******** Inference algorithm ********/
    std::cout << "Running inference algorithm..." << endl;
    time_t t_ini_abs = time(NULL);
    time_t t_ini;
    time_t t_end;
    bool stop = false;
    param.it = 0;
    int why = -1;
    double val_llh;
    int duration;
    while(!stop) {
        std::cout << " +Iteration " << param.it << "..." << endl;
        t_ini = time(NULL);    // Measure time elapsed

        // Run an inference step
        if(param.flag_ppca) {
            // TODO
            //my_infer_ppca::inference_step(data,hyper,param,pvar,semilla);
        } else {
            my_infer::inference_step(data,hyper,param,pvar,semilla);
        }

        // Write time elapsed
        t_end = time(NULL);    // Measure time elapsed
        my_output::write_telapsed(param.outdir,param.it,t_ini,t_end);
        duration = difftime(t_end,t_ini_abs);

        // Compute and write norm
        if(param.flag_userVec) {
            my_output::write_norm(param.outdir,param.it,my_infer::compute_avg_norm(pvar.rho),my_infer::compute_avg_norm(pvar.alpha),my_infer::compute_avg_norm(pvar.theta));
        } else {
            my_output::write_norm(param.outdir,param.it,my_infer::compute_avg_norm(pvar.rho),my_infer::compute_avg_norm(pvar.alpha));
        }

        // Compute llh + Check stop criteria
        if((param.it%param.rfreq)==0) {
            if(param.flag_ppca) {
                // TODO
                //val_llh = my_infer_ppca::compute_val_likelihood(true,duration,data,param,hyper,pvar);
                //my_infer_ppca::compute_test_performance(duration,data,hyper,param,pvar);
            } else {
                val_llh = my_infer::compute_val_likelihood(true,duration,data,param,hyper,pvar);
                my_infer::compute_test_performance(duration,data,hyper,param,pvar);
            }
            if(!param.noVal && (param.it>3*param.rfreq)) {
                /*
                std::cout << " (validation llh relative change (abs value): " << fabs((val_llh-param.prev_val_llh)/param.prev_val_llh) << ")" << endl;
                // If log likelihood increased, is not zero, and it increased less than 0.000001 of the previous value, set why to zero
                if(val_llh>=param.prev_val_llh && param.prev_val_llh!=0 && fabs((val_llh-param.prev_val_llh)/param.prev_val_llh)<param.valTolerance) {
                    stop = true;
                    why = 0;
                } else if(val_llh<param.prev_val_llh) {
                    // Count the number of times in a row that the likelihood decreased
                    (param.n_val_decr)++;
                } else if(val_llh>param.prev_val_llh) {
                    param.n_val_decr = 0;
                }
                if(param.n_val_decr>param.valConsecutive) {
                    stop = true;
                    why = 1;
                }
                */
            }
            param.prev_val_llh = val_llh;
        }
        if(param.it+1>=param.Niter) {
            stop = true;
            why = 2;
        }
        if(!stop) {
            (param.it)++;
            // Save results to file if desired
            if((param.it%param.saveCycle)==0) {
                my_output::write_all_matrices("_it"+std::to_string(param.it),data,hyper,param,pvar);
            }
        }
    }

    /******** Print output ********/
    // Compute final validation llh & Print max file
    if(param.flag_ppca) {
        // TODO
        //val_llh = my_infer_ppca::compute_val_likelihood(false,duration,data,param,hyper,pvar);
        //my_output::write_max_file(param,duration,val_llh,why);
        //my_infer_ppca::compute_test_performance(duration,data,hyper,param,pvar);
    } else {
        val_llh = my_infer::compute_val_likelihood(true,duration,data,param,hyper,pvar);
        my_output::write_max_file(param,duration,val_llh,why);
        my_infer::compute_test_performance(duration,data,hyper,param,pvar);
    }
    // Print output matrices
    my_output::write_all_matrices("",data,hyper,param,pvar);

    /******** Free memory ********/
    gsl_rng_free(semilla);

    /******** Return ********/
    return 0;
};

