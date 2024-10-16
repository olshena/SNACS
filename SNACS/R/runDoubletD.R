####################################################################
####################################################################
#' Run doubletD
#'
#' Assign singlet or doublet designation to each cell based on total depth and alternate depth values by applying doubletD method. The result is stored in the column "doubletD" in the "annCell" table of the SNACS object.
#'
#' @param snacsObj SNACSList object. Either snacsObj or depthTotal and depthAlt has to be provided
#' @param depthTotal Numeric matrix. Matrix having total depth values as numeric data. Rows depict SNPs and columns cells
#' @param depthAlt Numeric matrix. Matrix having total alternate values as numeric data. Rows depict SNPs and columns cells
#' @return A SNACSList object if input is a SNACSList; otherwise, a vector of singlet or doublet designation for each cell
#' @export
runDoubletD=function(snacsObj=NULL,depthTotal=NULL,depthAlt=NULL) {

    ###Set functions for doubletD
    getBetaMOM=function(x) {
        m_x=mean(x)
        minus_x=1- m_x
        s_x=stats::sd(x)
        y=((m_x * minus_x/ s_x^2) - 1)
        x_alpha=m_x * y
        x_beta=minus_x * y
        denom=x_alpha + x_beta
        ##  x_alpha=m_x * ((m_x * (1 - m_x) / s_x^2) - 1)
        ##  x_beta=(1 - m_x) * ((m_x * (1 - m_x) / s_x^2) - 1)

        ##    return(c(x_alpha / (x_alpha + x_beta), x_alpha + x_beta))
        return(c(x_alpha / denom, denom))
    }

    ###Sigma is global.  Theta is global.  Both are lists with two elements.  muts is a global vector
    ###Call to py_xz should be examamined.  As well cacall to px_z.
    ###py_xz is not a list but called like a list.  px_z is a list.  Both seem to have weird matching on variable names.
    ###The z argument is constant so it should stop being called
    ###The number of loops is very small, just mut*Sigma[[z]]*Theta[[z]]
    prv_z=function(cell,z,alpha_fp,alpha_fn,precision) {
        zminusone=z-1
        log_prob_sum=0
        
        ## --------------------------
        # Replaces function prv_y_bb to improve performance
        
        ## --------------------------

        for (mut in muts) {
            v=df_alt[cell, mut]
            #print(v)
            r=df_total[cell, mut] - v
            #print(r)
            prob_sum=0
            
            ## --------------------------
            # Replaces function prv_y_bb to improve performance
            
            n=r + v
            ## --------------------------

            for (x in Sigma[[z]]) {
                prob_sum_x=0
                for (y in Theta[[z]]) {
                    ## --------------------------
                    # Replaces function prv_y_bb to improve performance
                    yprime=alpha_fp + (1 - alpha_fp - alpha_fn) * y
                    if (yprime == 0) {
                      yprime=0.001
                    }
                    if (yprime == 1) {
                      yprime=0.999
                    }
                    alpha=precision * yprime
                    beta=precision - alpha

                    num=lgamma(n + 1) + lgamma(v + alpha) + lgamma(n - v + beta) + lgamma(alpha + beta)
                    den=lgamma(v + 1) + lgamma(n - v + 1) + lgamma(n + alpha + beta) + lgamma(alpha) + lgamma(beta)
                    prob_y_bb=exp(num - den)

                    ## --------------------------
                    
                    prob_sum_x=prob_sum_x + prob_y_bb * pyxz.matrix[which(round(y,3)==pyxz.matrix[,1] & x==pyxz.matrix[,2] & zminusone==pyxz.matrix[,3]),4]
                    ## prob_sum_x=prob_sum_x + prv_y_bb(r, v, y, alpha_fp, alpha_fn, precision) * pyxz.matrix[which(round(y,3)==pyxz.matrix[,1] & x==pyxz.matrix[,2] & zminusone==pyxz.matrix[,3]),4]
                    ##                prv_y_bb(r, v, y, alpha_fp, alpha_fn, precision) * py_xz[[paste(formatC(y, digits = 3, format = "f"),
                    ##                                  formatC(x, digits = 3, format = "f"), (z-1), sep = "_")]]
                    ##                prv_y_bb(r, v, y, alpha_fp, alpha_fn, precision) * py_xz[[paste(formatC(y, digits = 3, format = "f"),
                    ##                                  formatC(x, digits = 3, format = "f"), zminusone, sep = "_")]]
                }
                ##            prob_sum=prob_sum + as.numeric(pxz.matrix[which(mut==pxz.matrix[,1] & x==pxz.matrix[,2] & zminusone==pxz.matrix[,3]) ,4]) * prob_sum_x
                prob_sum=prob_sum + px_z[[paste(mut,x,zminusone,sep = "_")]] * prob_sum_x
                ##            prob_sum=prob_sum + px_z[[paste(mut,x,(z-1),sep = "_")]] * prob_sum_x
            }
            log_prob_sum=log(prob_sum) + log_prob_sum
        }
        return(log_prob_sum)
    }

    
    if (!is.null(snacsObj)) {
        df_total=snacsObj$depthTotal
        df_alt=snacsObj$depthAlt
    } else {
        df_total=depthTotal
        df_alt=depthAlt
    }
    if (is.null(df_total) | is.null(df_alt)) stop("Total depth and alternate depth matrices have to be provided")

    df_total=t(df_total)
    df_alt=t(df_alt)
    df_total=data.frame(df_total)
    df_alt=data.frame(df_alt)

    ## ####Index input data by cell
    ## rownames(df_total)=df_total[,1]
    ## df_total=df_total[, -1]

    ## rownames(df_alt)=df_alt[,1]
    ## df_alt=df_alt[, -1]

    cells=rownames(df_total)
    muts=colnames(df_total)

    ###Calculate VAF, drops VAFs with divide by zero error and turns matrix into a vector
    vaf_values=as.matrix(df_alt/df_total)
    #vaf_values=df_alt/df_total
    vaf_values=vaf_values[!is.na(vaf_values) & is.finite(vaf_values)]

    ###Set parameters for doubletD
    delta=0.2
    beta=0.05
    Sigma=list(c(0, 1/2, 1), c(0, 1/4, 1/2, 3/4, 1))
    Theta=list(c(0, 1/2, 1), c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1))
    threshold=log((1 - delta) / delta)

    ###Calc Alpha false neg and false pos, Precision
    mean1_prec1=getBetaMOM(vaf_values[which(vaf_values <= 0.15)])
    alpha_fp=mean1_prec1[1]

    mean2_prec2=getBetaMOM(vaf_values[which(vaf_values >= 0.85)])
    alpha_fn=1 - mean2_prec2[1]

    mean3_prec3=getBetaMOM(vaf_values[which((vaf_values > 0.15) & (vaf_values < 0.85))])
    precision=stats::median(c(mean1_prec1[2], mean2_prec2[2], mean3_prec3[2]))

    ###Set px_z
    #####Initialize px_z dictionary-like structure as a list
    ###Seems to be a vector structured as a list
    px_z=c()
    ##px_z=list()

    for (z in seq_along(Sigma)) {
      for (m in muts) {
        for (sigma in Sigma[[z]]) {
          px_z_key=paste(m, sigma, (z-1), sep = "_")
          px_z[[px_z_key]]=0
        }
      }
    }

    for (m in muts) {
        reads=data.frame(df_alt[[m]])
        rownames(reads)=rownames(df_alt)
        colnames(reads)=m
        total=data.frame(df_total[[m]])
        rownames(total)=rownames(df_total)
        colnames(total)=m
        vaf=data.frame(reads/total)

        loh_cells=vaf[which(vaf[, m] > 0.85), , drop = FALSE]
        loh_cells=stats::na.omit(loh_cells)

        wt_cells=vaf[which(vaf[, m] < 0.15), , drop = FALSE]
        wt_cells=stats::na.omit(wt_cells)

        het_cells=vaf[which(vaf[, m] < 0.85 & vaf[, m] > 0.15), , drop = FALSE]
        het_cells=stats::na.omit(het_cells)

        non_na_cells =  dim(loh_cells)[1] + dim(wt_cells)[1] + dim(het_cells)[1]
        est_wt_rate = dim(wt_cells)[1]/non_na_cells
        est_loh_rate = dim(loh_cells)[1]/non_na_cells
        est_het_rate = dim(het_cells)[1]/non_na_cells

        px_z[paste(m, 0, 0, sep = "_")] = est_wt_rate
        px_z[paste(m, 0.5, 0, sep = "_")] = est_het_rate
        px_z[paste(m, 1, 0, sep = "_")] = est_loh_rate
        }

    for (m in muts) {
        for (a in Sigma[[1]]) {
            for (b in Sigma[[1]]) {
                cc=(a + b) / 2
                if (cc %in% Sigma[[2]]) {
                    px_z[[paste(m, cc, 1, sep = "_")]] = px_z[[paste(m, cc, 1, sep = "_")]] +
                    px_z[[paste(m, a, 0, sep = "_")]] * px_z[[paste(m, b, 0, sep = "_")]]
                }
            }
        }
    }

    ##pxz.matrix=cbind(matrix(unlist(strsplit(names(px_z),split="_")),ncol=3,byrow=TRUE),as.numeric(px_z))

    ####Set norm_const
    norm_const=list()

    for (m in muts) {
        sum_value=0


      for (a in Sigma[[1]]) {
        for (b in Sigma[[1]]) {
          value_name_a=paste(m, a, 0, sep = "_")
          value_name_b=paste(m, b, 0, sep = "_")
          sum_value=sum_value + px_z[[value_name_a]] * px_z[[value_name_b]]
        }
      }
      norm_const[m]=sum_value
    }

    for (m in muts) {
      for (cc in Sigma[[1]]) {
        px_z[[paste(m, cc, 1, sep = "_")]]=px_z[[paste(m, cc, 1, sep = "_")]] / norm_const[[m]]
      }
    }

    ###Set py_xy
    ##### Initialize py_xy dictionary-like structure as a list
    py_xz=c()
    for (z in c(0, 1)) {
      for (theta in Theta[[z + 1]]) {
        for (sigma in Sigma[[z + 1]]) {
          key=paste0(formatC(theta, digits = 3, format = "f"),
                        "_", formatC(sigma, digits = 3, format = "f"),
                        "_", z)
          py_xz[key]=0
        }
      }
    }

    #####Calculate py_xz for all combinations and add to list
    py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.000, digits = 3, format = "f"), 0, sep = "_")]] = 1 - beta**2

    py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 0, sep = "_")]] = beta * (1 - beta)

    py_xz[[paste(formatC(0.500, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 0, sep = "_")]] = (1-beta)**2

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 0, sep = "_")]] = beta * (1 - beta)

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(1.000, digits = 3, format = "f"), 0, sep = "_")]] = 1 - beta**2

    py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.000, digits = 3, format = "f"), 1, sep = "_")]] = 1 -beta**4

    py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]] = beta * (1 - beta)**3 +
                 3 * beta**2 * (1 - beta)**2 + 3 * beta**3 * (1 - beta)

    py_xz[[paste(formatC(0.250, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]] = (1 - beta)**4

    py_xz[[paste(formatC(0.333, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]] = 3*beta*(1 - beta)**3

    py_xz[[paste(formatC(0.500, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]] = 3*beta**2 * (1 - beta)**2

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]] = beta**3 * (1 - beta)

    py_xz[[paste(formatC(0.333, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 1, sep = "_")]] = 2 * beta * (1 - beta)**3

    py_xz[[paste(formatC(0.500, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 1, sep = "_")]] = (1-beta)**4 +
                                                               4 * beta**2 * (1 - beta)**2
    py_xz[[paste(formatC(0.667, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 1, sep = "_")]] = 2 * beta * (1 - beta)**3

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(0.500, digits = 3, format = "f"), 1, sep = "_")]] = beta**2 * (1 - beta)**2 + 2 * beta**3 * (1 - beta)

    py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.750, digits = 3, format = "f"), 1, sep = "_")]] = py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]]

    py_xz[[paste(formatC(0.500, digits = 3, format = "f"),
                 formatC(0.750, digits = 3, format = "f"), 1, sep = "_")]] = py_xz[[paste(formatC(0.500, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]]

    py_xz[[paste(formatC(0.667, digits = 3, format = "f"),
                 formatC(0.750, digits = 3, format = "f"), 1, sep = "_")]] = py_xz[[paste(formatC(0.333, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]]

    py_xz[[paste(formatC(0.750, digits = 3, format = "f"),
                 formatC(0.750, digits = 3, format = "f"), 1, sep = "_")]] = py_xz[[paste(formatC(0.250, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]]

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(0.750, digits = 3, format = "f"), 1, sep = "_")]] = py_xz[[paste(formatC(0.000, digits = 3, format = "f"),
                 formatC(0.250, digits = 3, format = "f"), 1, sep = "_")]]

    py_xz[[paste(formatC(1.000, digits = 3, format = "f"),
                 formatC(1.000, digits = 3, format = "f"), 1, sep = "_")]] = 1 - beta**4

    pyxz.matrix=cbind(matrix(as.numeric(unlist(strsplit(names(py_xz),"_"))),ncol=3,byrow=TRUE),as.numeric(py_xz))

    ###Calculate Log Prob-- this takes a really long time!
    logprobs_1=rep(NA,length(cells))
    for (i in 1:length(cells)) {
      logprobs_1[i]=prv_z(cell=cells[i],z=1,alpha_fp=alpha_fp,alpha_fn=alpha_fn,precision=precision)
    }

    logprobs_2=rep(NA,length(cells))
    for (i in 1:length(cells)) {
      logprobs_2[i]=prv_z(cell=cells[i],z=2,alpha_fp=alpha_fp,alpha_fn=alpha_fn,precision=precision)
    }

    #rownames(logprobs_1)=rownames(df_total)
    #rownames(logprobs_2)=rownames(df_total)

    #logprobs=cbind(logprobs_1, logprobs_2)

    ## Calculate the final call as doublet vs singlet and output to doublet_result vector
    doublet_result=ifelse(logprobs_2 - logprobs_1 > threshold, "Doublet", "Singlet")
    names(doublet_result)=rownames(df_total)
    
    if (!is.null(snacsObj)) {
        snacsObj$annCell$doubletD=doublet_result
        invisible(snacsObj)
    } else {
        invisible(doublet_result)
    }
}

####################################################################
####################################################################
#' Add doubletD calls to SNACS calls
#'
#' First run doubletD, which assigns singlet or doublet designation to each cell based on total depth and alternate depth values. Then designate those cells as a doublet that are called so by doubletD but are called a singlet by SNACS Round 2. The result is stored in the column "snacsPlusDoubletD" in the "annCell" table of the SNACS object.
#'
#' @param snacsObj SNACSList object
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
runSNACSplusDoubletD=function(snacsObj,verbose=FALSE) {
    if (is.na(match("snacsRnd1",names(snacsObj$annCell)))) stop("Run makeHashCall() before calling runSNACSplusDoubletD()\n")
    
    snacsObj=runDoubletD(snacsObj)
    
    if ("snacsRnd2"%in%names(snacsObj$annCell)) {
        snacsObj$annCell$snacsPlusDoubletD=snacsObj$annCell$snacsRnd2
    } else {
        snacsObj$annCell$snacsPlusDoubletD=snacsObj$annCell$snacsRnd1
    }
    snacsObj$annCell$snacsPlusDoubletD[which(snacsObj$annCell$snacsPlusDoubletD%in%snacsObj$annHash$hashNames & snacsObj$annCell$doubletD=="Doublet")]="Doublet"

    if (verbose) cat('"doubletD" and "snacsPlusDoubletD" columns added to "annCell" table in SNACS object\n\n',sep="")
    
    invisible(snacsObj)
}

###########################################################
###########################################################
