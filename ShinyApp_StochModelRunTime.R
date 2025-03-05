## R shiny app to test computer run time of the 'time-to-next-event' Algorithm (aka Gillaspie Algorithm):
## Model: Stochastic EBOLA model per Legrand et al. 2007
## 4/8/18, by Wan Yang

library(shiny)
library(deSolve)



ui <- fluidPage(h2('Computing Time'),
                
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'Nwk',h6('Simulation Length (Week)'),value=9,min=9,max=17,step=1),
                               h6('Ebola Model per the Gillaspie Algorithm (Legrand et al. 2007'),
                               h6('WAIT...WAIT...WAIT...'),
                               h6('IT MAY TAKE A WHILE...')
                               ),
                              
                              mainPanel(
                                plotOutput(outputId = 'plots',width = '100%', height = "550px")
                                )
                              )
                
                )

server <- function(input, output){
  
  output$plots=renderPlot({
    
    num_wk=input$Nwk  # number of week to run
    
    alpha=1; # incubation period: 7 days = 1 week
    gamma.h=7/5; # from onset to hopspitalization: 5 days
    gamma.d=7/9.6; # from onset to death: 9.6 days
    gamma.i=7/10; # from onset to end of infectiousness for survivors: 10 days
    gamma.f=7/2; # from death to traditional burial 2 days
    gamma.ih=7/(10-5); # from hospitalization to end of infectiousness for survivors
    gamma.dh=7/(9.6-5); # from hospitalization to death
    theta1=.67; # proportion infectious in the hospital
    delta1=.8; # CFR for unhospitalized
    delta2=.8; # CFR for hospitalize
    ## reduced by a factor of p as we are running it deterministically
    beta.I=.588; # transmission rate in the community
    beta.H=.794; # transmission rate in the hospital
    beta.F=7.653; # transmission rate at funerals; 
    
    # INTIAL CONDITIONS
    repeat {
      N=2e5; E0=H0=FF0=R0=0; I0=3; S0=N-I0; 
      tm=0; # in week
      S=S0; E=E0; I=I0; H=H0; FF=FF0; R=R0; cumI=I0;
      res=c(tm,NA,S,E,I,H,FF,R,cumI); 
      
      tm1=proc.time()["user.self"] # record time when it begins
      
      while((I>0 | E>0 | H>0 | FF>0) & S>0 & tm<num_wk){ # simulate epid for the first 20 week without ctrl
        # step 1: compute the transition rates
        rate.se=S/N*(beta.I*I+beta.H*H+beta.F*FF);  # S->E
        rate.ei=alpha*E;  # E->I
        rate.ih=gamma.h*theta1*I; # I->H
        rate.hf=gamma.dh*delta2*H;  # H->F
        rate.fr=gamma.f*FF;  # F->R
        rate.ir=gamma.i*(1-theta1)*(1-delta1)*I;  # I->R
        rate.if=delta1*(1-theta1)*gamma.d*I; # I->F
        rate.hr=gamma.ih*(1-delta2)*H; # H->R
        rate.tot=rate.se+rate.ei+rate.ih+rate.hf+rate.fr+rate.ir+rate.if+rate.hr
        
        # step 2: draw a random number, u1, betw 0 and 1 and 
        # calcualte the time after which the next transition occurs
        u1=runif(1); # draw a random number 
        tau_i=-log(u1)/rate.tot
        
        # step 3: compute the probability that each type of transition will occur based on the rates
        # use this to calculate the range in which a number drawn at random 
        # must lie for a given transition to occur
        p.se=rate.se/rate.tot;  # prob the next event is S->E
        p.ei=rate.ei/rate.tot;
        p.ih=rate.ih/rate.tot;
        p.hf=rate.hf/rate.tot;
        p.fr=rate.fr/rate.tot;
        p.ir=rate.ir/rate.tot;
        p.if=rate.if/rate.tot;
        p.hr=rate.hr/rate.tot;
        
        # step 4: draw a random number, u2, to determine the transition event which occurs next.
        u2=runif(1); # draw a random number 
        
        # step 5: use the result from step 4 to update the number of state varibles
        if(u2<p.se) { # u2 lies in [0,p.se), the next event is S->E
          S=S-1; E=E+1;
        } else if (u2<p.se+p.ei) { # u2 lies in [p.se,p.se+p.ei), the next event is E->I
          E=E-1; I=I+1;
          cumI=cumI+1; # record cumulative incidence
        } else if (u2<p.se+p.ei+p.ih) { # u2 lies in [p.se+p.ei,p.se+p.ei+p.ih), the next event is I->H
          I=I-1; H=H+1;
        } else if (u2<p.se+p.ei+p.ih+p.hf) { # u2 lies in [p.se+p.ei+p.ih,p.se+p.ei+p.ih+p.hf), the next event is H->F
          H=H-1; FF=FF+1;
        } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr){ # FF->R
          FF=FF-1; R=R+1;
        } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir){ # I->R
          I=I-1; R=R+1;
        } else if (u2<p.se+p.ei+p.ih+p.hf+p.fr+p.ir+p.if) { # I->FF
          I=I-1; FF=FF+1;
        } else { # H->R
          H=H-1; R=R+1;
        }
        tm=tm+tau_i; # update time
        # save the results
        res=rbind(res,c(tm,tau_i,S,E,I,H,FF,R,cumI))
      } 
      colnames(res)=c('time','tm_step','S','E','I','H','FF','R','cumI')
      
      tm_span=proc.time()["user.self"]-tm1  # computing time
      
      # to make sure the epidemic does take off
      if(tail(res[,'cumI'],1)>exp(1.25+.42*num_wk)*.8) break;  
    }
    
    
    # compile number of newI, I, H, FF, by week
    wkly.res=matrix(NA,num_wk,8); colnames(wkly.res)=c('time','S','E','I','H','FF','R','cumI')
    num_events=NULL;
    for(wk in 1:num_wk){
      idx=tail(which(res[,'time']<wk),1);
      # idx also records the number of events happened from time 0 to the current week
      num_events=c(num_events,ifelse(wk==1,idx,idx-sum(num_events)));
      wkly.res[wk,]=res[idx,c('time','S','E','I','H','FF','R','cumI')];
    }
    
    
    # plot results
    par(mfrow=c(2,1),mar=c(3,3,1,1),cex=1.2,mgp=c(1.6,.5,0),cex.axis=.9)
    plot(1:num_wk,wkly.res[,'cumI'],xlab='Week',ylab='Cumulative Incidence',lwd=2,pch='x',
         main=paste('Run time:',format(tm_span,digits = 2),'seconds'),cex.main=1.5)
    matplot(1:num_wk,wkly.res[,c('E','I','H','FF')],xlab='Week',ylab='Number of Ebola-related People',col=rainbow(4),type='l',lty=1,lwd=2)
    legend('topleft',legend = c('Exposed','Infectious in the community','Infectious in the hospitical','Deceased, not yet buried'),
           col=rainbow(4),lty=1,lwd=2,bty='n')
  })
  
}

shinyApp(ui=ui, server = server)