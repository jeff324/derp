
Model = function(lp,pars) {

     model_ind = R6::R6Class(
          'Model.Individual',
          public = list(
               lp = NULL,
               theta = NULL,
               pars = NULL,
               initialize = function(pars,lp) {
                    self$lp = lp
                    self$theta = names(pars)
                    self$pars = pars
               },
               initializer = function(par_name) {
                    lb = self$pars[[par_name]]$init[1]
                    ub = self$pars[[par_name]]$init[2]
                    rv = runif(1,lb,ub)
                    return(rv)
               }
          )
     )

     model_hier = R6::R6Class(
          'Model.Hierarchical',
          inherit = model_ind,
          public = list(
               lp = NULL,
               phi = NULL,
               blocks = NULL,
               levels = NULL,
               initialize = function(pars,lp) {
                    self$lp = lp
                    self$levels = private$get_levels(pars)
                    self$theta = private$get_pars(pars,level=1)
                    self$phi = private$get_pars(pars,level=2)
                    self$pars = pars
                    self$blocks = private$get_blocks(pars)
               }
          ),
          private = list(

               get_levels = function(pars)
               {
                    level_idx = sapply(pars,function(x) x$block)
                    level_idx = sapply(level_idx,function(x)
                    {
                         if (is.null(x))
                         {
                              x = 1
                         } else if (x == 0)
                         {
                              x = 1
                         } else if (x != 0)
                         {
                              x = 2
                         }
                         return(x)
                    })
                    return(level_idx)
               },

               get_blocks = function(pars)
               {
                    pars = pars[which(self$levels==2)]
                    block = sapply(1:length(pars),function(i) pars[[i]]$block)
                    block_list = list()
                    unique_blocks = unique(block)
                    idx = 1:length(block)
                    block_list = lapply(unique_blocks,function(x) idx[block == x])
                    return(block_list)
               },

               get_pars = function(pars,level)
               {
                    names(pars)[which(self$levels == level)]
               }

          )
     )


     if (any(!is.null(unlist(sapply(pars,function(x) x$block))))) {
          model = model_hier$new(pars,lp)
     } else {
          model = model_ind$new(pars,lp)

     }

     return(model)

}



