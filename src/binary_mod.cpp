#ifndef MB_BASE_PLUGIN
#define MB_BASE_PLUGIN <optionals.h>
#endif

#include <binary.h>


RCPP_MODULE(binary){  
  Rcpp::class_< MB_Base<binary> >("binary")
    
    .constructor<const List>()
    
    .method( "get.f", & MB_Base<binary>::get_f)
    .method( "get.df", & MB_Base<binary>::get_df)
    .method( "get.fdf", & MB_Base<binary>::get_fdf)
    .method( "get.fdfh", & MB_Base<binary>::get_fdfh)
    .method( "get.tape.stats", & MB_Base<binary>::get_tape_stats)
    .method( "get.hessian", & MB_Base<binary>::get_hessian)
    .method( "get.hessian.sparse", & MB_Base<binary>::get_hessian_sparse)
    .method( "record.tape", & MB_Base<binary>::record_tape)
    .method( "init.hessian", & MB_Base<binary>::hessian_init)
    .method( "init.sparse.hessian", & MB_Base<binary>::hessian_init)
    .method( "get.f.direct", & MB_Base<binary>::get_f_direct) // plugin
    .method( "get.LL", & MB_Base<binary>::get_LL) // plugin
    ;
}

