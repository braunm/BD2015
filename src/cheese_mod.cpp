#ifndef MB_BASE_PLUGIN
#define MB_BASE_PLUGIN <optionals.h>
#endif

#include <cheese.h>


RCPP_MODULE(cheese){  
  Rcpp::class_< MB_Base<cheese> >("cheese")
    
    .constructor<const List>()
    
    .method( "get.f", & MB_Base<cheese>::get_f)
    .method( "get.df", & MB_Base<cheese>::get_df)
    .method( "get.fdf", & MB_Base<cheese>::get_fdf)
    .method( "get.fdfh", & MB_Base<cheese>::get_fdfh)
    .method( "get.tape.stats", & MB_Base<cheese>::get_tape_stats)
    .method( "get.hessian", & MB_Base<cheese>::get_hessian)
    .method( "get.hessian.sparse", & MB_Base<cheese>::get_hessian_sparse)
    .method( "record.tape", & MB_Base<cheese>::record_tape)
    .method( "init.hessian", & MB_Base<cheese>::hessian_init)
    .method( "init.sparse.hessian", & MB_Base<cheese>::hessian_init)
    .method( "get.f.direct", & MB_Base<cheese>::get_f_direct) // plugin
    .method( "get.LL", & MB_Base<cheese>::get_LL) // plugin
    ;
}

