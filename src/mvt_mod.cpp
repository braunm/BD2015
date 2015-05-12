#ifndef MB_BASE_PLUGIN
#define MB_BASE_PLUGIN <optionals.h>
#endif

#include <mvt.h>


RCPP_MODULE(mvt){  
  Rcpp::class_< MB_Base<mvt> >("mvt")
    
    .constructor<const List>()
    
    .method( "get.f", & MB_Base<mvt>::get_f)
    .method( "get.df", & MB_Base<mvt>::get_df)
    .method( "get.fdf", & MB_Base<mvt>::get_fdf)
    .method( "get.fdfh", & MB_Base<mvt>::get_fdfh)
    .method( "get.tape.stats", & MB_Base<mvt>::get_tape_stats)
    .method( "get.hessian", & MB_Base<mvt>::get_hessian)
    .method( "get.hessian.sparse", & MB_Base<mvt>::get_hessian_sparse)
    .method( "record.tape", & MB_Base<mvt>::record_tape)
    .method( "init.hessian", & MB_Base<mvt>::hessian_init)
    .method( "init.sparse.hessian", & MB_Base<mvt>::hessian_init)
    .method( "get.f.direct", & MB_Base<mvt>::get_f_direct) // plugin
    ;
}

