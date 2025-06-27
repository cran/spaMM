if (spaMM.getOption("example_maxtime")>4) {
  cat(cli::col_yellow("\ntest-mm:"))
  data("Gryphon")
  noPAIR <- Gryphon_df
  Gryphon_df$PAIR <- 1 # should be 1 
  
  (byPAIR <- fitme(BWT ~ 1+ corrMatrix(PAIR | PAIRfn(ID,mother)), 
        data=Gryphon_df, corrMatrix=Gryphon_A, control.HLfit=list(algebra="spprec"))) # -2023.216
  
  (bymm <- fitme(BWT ~ 1+ corrMatrix(mm(ID,mother) | mmfn(ID,mother)), 
        data=noPAIR, corrMatrix=Gryphon_A, control.HLfit=list(algebra="spprec"))) # -2023.216
  
  testthat::expect_true(diff(range(logLik(byPAIR),logLik(bymm),-2023.215998842957))<1e-12) 
  
  
  if (FALSE) {
    
    (byPAIRsx <- fitme(BWT ~ 1+ corrMatrix(PAIR*sex | PAIRfn(ID,mother)), 
                    # fixed=list(ranCoefs=list("1"=c(1,0,0,0,1,0,0,1,0,1)),phi=4),
                    data=Gryphon_df, corrMatrix=Gryphon_A, control.HLfit=list(algebra="spprec"))) # -1977.144
    
    (bymmsx <- fitme(BWT ~ 1+ corrMatrix(mm(ID,mother)*sex | mmfn(ID,mother)), 
                  # fixed=list(ranCoefs=list("1"=c(1,0,0,0,1,0,0,1,0,1)),phi=4),
                  data=noPAIR, corrMatrix=Gryphon_A, control.HLfit=list(algebra="spprec"))) # -1977.144
    
    
    testthat::expect_true(diff(range(logLik(byPAIRsx),logLik(bymmsx),-1977.143641255906))<1e-12) 
    
    
  }
  
  if (FALSE) { # 0 + not the same.  Should check factor.[asgn[colit]] in such cases.
    
    # reversed results with new PAIRfn() ie mmfn()
    
    fitme(BWT ~ 1+ corrMatrix(0+PAIR | PAIRfn(ID,mother)),
          data=Gryphon_df, corrMatrix=Gryphon_A) # old PAIRfn: stop()
    
    fitme(BWT ~ 1+ corrMatrix(0+mm(ID,mother) | PAIRfn(ID,mother)),
          data=Gryphon_df, corrMatrix=Gryphon_A) # old PAIRfn: different likelihood 
    
    
    
    
    # fitme(BWT ~ 1+ corrMatrix(0+mm(ID,mother)*sex | PAIRfn(ID,mother)),
    #      data=Gryphon_df, corrMatrix=Gryphon_A)
    
    
    # fitme(BWT ~ 1+ corrMatrix(0+mm(ID,mother) | PAIRfn(ID,mother)),
    #      data=noPAIR, corrMatrix=Gryphon_A)
    
    # fitme(BWT ~ 1+ corrMatrix(0+mm(ID,mother)*sex | PAIRfn(ID,mother)),
    #       data=noPAIR, corrMatrix=Gryphon_A)
    
    
  }
}
