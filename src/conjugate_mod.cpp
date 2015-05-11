#ifndef MB_BASE_PLUGIN
#define MB_BASE_PLUGIN <optionals.h>
#endif

#include <conjugate.h>


RCPP_MODULE(conjugate){  
  Rcpp::class_< MB_Base<conjugate> >("conjugate")
    
    .constructor<const List>()
    
    .method( "get.f", & MB_Base<conjugate>::get_f)
    .method( "get.df", & MB_Base<conjugate>::get_df)
    .method( "get.fdf", & MB_Base<conjugate>::get_fdf)
    .method( "get.fdfh", & MB_Base<conjugate>::get_fdfh)
    .method( "get.tape.stats", & MB_Base<conjugate>::get_tape_stats)
    .method( "get.hessian", & MB_Base<conjugate>::get_hessian)
    .method( "get.hessian.sparse", & MB_Base<conjugate>::get_hessian_sparse)
    .method( "record.tape", & MB_Base<conjugate>::record_tape)
    .method( "init.hessian", & MB_Base<conjugate>::hessian_init)
    .method( "init.sparse.hessian", & MB_Base<conjugate>::hessian_init)
    .method( "get.f.direct", & MB_Base<conjugate>::get_f_direct) // plugin
    .method( "get.LL", & MB_Base<conjugate>::get_LL) // plugin
    //   .method( "get.LLi", & MB_Base<conjugate>::get_LLi) // plugin
    //   .method( "get.prior.i", & MB_Base<conjugate>::get_prior_i) // plugin
    ;
}

