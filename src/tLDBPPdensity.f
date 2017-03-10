c=======================================================================                      
      subroutine tldbppdensity(y,x,nrec,p,
     &                         npred,ngrid,grid,xpred,
     &                         maxn,nu,alpha,lambda,tau1,tau2,psiinv,
     &                         s0invm,s0inv,
     &                         kk,gp,beta,theta,mub,sb,
     &                         mcmc,nsave,slice,
     &                         acrate,thetasave,randsave,
     &                         fmean,flow,fupp,meanfpm,meanfpl,meanfph, 
     &                         cpo,seed,
     &                         iflag,sbinv,workm1,workv1,workv2,
     &                         workmh1,workmh2,betal,betar,beta1,
     &                         thetal,thetar,theta1,workdpw,weight,
     &                         fw,fw2,fs,fm,worksam,v)

c=======================================================================
c     # 58 arguments
c
c     Subroutine `tldbppdensity' to run a Markov chain for a 
c     singles atoms Linear Dependent Bernstein-Dirichlet Process  
c     prior for bounded conditional density estimation.
c
c     Copyright: Andrés F. Barrientos and Alejandro Jara, 2010.
c
c     Version 1.0: 
c
c     Last modification: 07-10-2012.
c     
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or (at
c     your option) any later version.
c
c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c     The author's contact information:
c
c      Felipe Barrientos
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Email: afbarrie@uc.cl
c
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
c      Fax  : +56-2-3547729  Email: atjara@uc.cl
c
c---- Data -------------------------------------------------------------
c 
c        nrec        :  integer giving the number of data points. 
c        p           :  integer giving the number of predictors.
c        y           :  real vector giving the responses, y(nrec). 
c        x           :  real matrix giving the design matrix, x(nrec,p).
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        cband       :  integer value indicating whether the 
c                       credible bands need to be computed or not.
c        ngrid       :  integer giving the number of grid points where 
c                       the density estimates are evaluated.. 
c        npred       :  integer giving the number of predictions.
c        grid        :  real vector giving the grid, grid(ngrid). 
c        tband       :  integer indicating the type of credible 
c                       band that need to be computed.
c        xpred       :  real matrix giving the design matrix of the 
c                       predictions, zpred(npred,p).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        alpha       :  real vector giving the parameters of the beta 
c                       distribution - the DP base measure -, alpha(2).
c        lambda      :  real giving the parameter of the truncated
c                       Poisson prior for the degree of the Bernstein
c                       polynomial.
c        maxn        :  integer giving the truncation of the DP.
c        tau1,tau2   :  real giving the values of hyperparameters of the
c                       Gamma prior for the g-prior (only internal use).
c        s0invm      :  real vector giving the product of the prior
c                       precision and the prior mean of the normal
c                       normal prior for mub, s0ivnm(p).
c        s0inv       :  real vector giving the inverse of the covariance
c                       matrix of the normal prior for mub, 
c                       s0inv(p,p).
c        nu          :  integer giving the degres of freedom parameter
c                       of the inverted-Wishart prior for sb.
c        psiinv      :  real matrix giving the scale matrix of the
c                       inverted-Wishart prior for sd, psiinv(p,p).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        kk          :  integer giving the degree of the Bernstein 
c                       polynomial.
c        gp          :  real giving the value of the G-prior 
c                       parameter (only for internal use).
c        beta        :  real matrix giving the value of the regression  
c                       coefficients, beta(maxn,p).
c        mub         :  real vector giving the mean of the normal 
c                       centering distribution, mub(p).
c        sb          :  real matrix giving the variance of the normal
c                       centering distribution, sb(p,p).
c        theta       :  real vector giving the atoms parameters,    
c                       theta(maxn).
c
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        nsave       :  integer giving the number of scans to be saved.
c        betaw       :  real giving the Slice sampling parameter for
c                       the regression coefficients.
c        gw          :  real giving the Slice sampling parameter for
c                       the g-prior parameter.
c        thetaw      :  real giving the Slice sampling parameter for
c                       the atoms parameters.
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the average number of
c                       evaluation of the posterior in the Slice
c                       sampling step and the acceptance rate 
c                       of the MH step for the degree of the 
c                       polynomial, acrate(2).
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the regression coeff and stick-breaking
c                       parameters, randsave(nsave,p*(maxn-1)+maxn).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,2+p+p*(p+1)/2)
c        fmean       :  real matrix giving the posterior mean of the 
c                       density, fmean(npred,ngrid).
c        flow        :  real matrix giving the lower limit of the 
c                       HPD of the density, flow(npred,ngrid).
c        fupp        :  real matrix giving the upper limit of the  
c                       HPD of the density, fupp(npred,ngrid).
c        meanfpm     :  real vector giving the posterior mean of the 
c                       mean function, meanfpm(npred).
c        meanfpl     :  real vector giving the lower limit of the 
c                       HPD of the mean function meanfpl(npred).
c        meanfph     :  real vector giving the upper limit of the  
c                       HPD of the mean function, meanfph(npred).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        betal       :  real matrix used for the slice-sampling of the
c                       regression coefficients, betal(maxn,p).
c        betar       :  real matrix used for the slice-sampling of the
c                       regression coefficients, betar(maxn,p).
c        beta1       :  real matrix used for the slice-sampling of the
c                       regression coefficients, beta1(maxn,p).
c        fs          :  real vector used to evaluate the conditional
c                       densities, fs(ngrid).
c        fm          :  real vector used to evaluate the conditional
c                       means, fs(npred).
c        fw          :  real vector used to compute cpos, 
c                       fw(nrec).
c        fw2         :  real matrix used to evaluate the predictions,
c                       fw2(npred,ngrid).
c        iflag       :  integer vector used to evaluate reg. coeff.,
c                       iflag(p).
c        sbinv       :  real matrix used to keep the inverse of sb,
c                       sbinv(p,p).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        thetal      :  real vector used for the slice sampling
c                       of atoms parameters, thetal(maxn).
c        thetar      :  real vector used for the slice sampling
c                       of atoms parameters, thetar(maxn).
c        theta1      :  real vector used for the slice sampling
c                       of atoms parameters, theta1(maxn).
c        weight      :  real vector used to compute the DP weights,
c                       weight(maxn).
c        workdpw     :  real vector used to compute the DP weights,
c                       workdpw(maxn+1).                       
c        workm1      :  real matrix used to update centering param.,
c                       workm1(p,p).
c        workmh1     :  real vector used to update centering param.,
c                       workmh1(p*(p+1)/2).
c        workmh2     :  real vector used to update centering param.,
c                       workmh2(p*(p+1)/2).
c        workv1      :  real vector used to update centering param.,
c                       workv1(p).
c        workv2      :  real vector used to update centering param.,
c                       workv2(p).
c        worksam     :  real vector used to comput HPD bands,
c                       worksam(nsave).
c        v           :  real vector giving the stick-breaking
c                       parameters, v(maxm).
c=======================================================================                  

      implicit none 

c+++++Data
      integer nrec,p
      real(kind=8) x(nrec,p)
      real(kind=8) y(nrec)

c+++++Predictions
      integer cband,npred,ngrid,tband
      real(kind=8)  grid(ngrid),xpred(npred,p)

c+++++Prior 
      integer maxn,nu
      real(kind=8) alpha(2)
      real(kind=8) lambda
      real(kind=8) tau1,tau2
      real(kind=8) psiinv(p,p)
      real(kind=8) s0invm(p)
      real(kind=8) s0inv(p,p)
      
c+++++Current values of the parameters
      integer kk
      real(kind=8) gp
      real(kind=8) beta(maxn-1,p)
      real(kind=8) theta(maxn)
      real(kind=8) mub(p)
      real(kind=8) sb(p,p)
      real(kind=8) v(maxn)


c+++++MCMC parameters
      integer mcmc(5),nburn,nskip,nsave,ndisplay
      real(kind=8) slice(3),betaw,gw,thetaw

c+++++Output
      real(kind=8) acrate(2)
      real(kind=8) thetasave(nsave,1+p+p*(p+1)/2+1)
      real(kind=8) randsave(nsave,p*(maxn-1)+maxn)
      real(kind=8) fmean(npred,ngrid)
      real(kind=8) flow(npred,ngrid)
      real(kind=8) fupp(npred,ngrid)
      real(kind=8) meanfpm(npred)
      real(kind=8) meanfpl(npred)
      real(kind=8) meanfph(npred)
      real(kind=8) cpo(nrec,2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ Inversion
      integer iflag(p)
      real(kind=8) sbinv(p,p)
      real(kind=8) workm1(p,p)
      real(kind=8) workv1(p)
      real(kind=8) workv2(p)
      real(kind=8) workmh1(p*(p+1)/2)
      real(kind=8) workmh2(p*(p+1)/2)

c++++ Slice sampler
      real(kind=8) betal(maxn-1,p)
      real(kind=8) betar(maxn-1,p)
      real(kind=8) beta1(maxn-1,p)
      real(kind=8) thetal(maxn)
      real(kind=8) thetar(maxn)
      real(kind=8) theta1(maxn) 

c++++ DP
      real(kind=8) workdpw(maxn+1)
      real(kind=8) weight(maxn)

c++++ CPO
      real(kind=8) fw(nrec)

c++++ predictions
      real(kind=8) fw2(npred,ngrid)

c++++ credible bands
      real(kind=8) fs(ngrid)
      real(kind=8) fm(npred) 
      real(kind=8) worksam(nsave)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer count,i,i1,ii,j,j1,k
      integer sprint
      real(kind=8) tmp1,tmp2

c+++++Distributions
      real(kind=8) dbet,dpoiss
      
c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++Slice sampler
      integer kk1
      real(kind=8) evaly,evalf,aux
      real(kind=8) gl
      real(kind=8) gr
      real(kind=8) g1

c+++++CPU time
      real(kind=8) sec00,sec0,sec1,sec

c+++++RNG and distributions
      real runif

c++++ opening files

      open(unit=1,file='dppackage1.out',status='unknown',
     &     form='unformatted')
      open(unit=2,file='dppackage2.out',status='unknown',
     &     form='unformatted')


c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      cband=mcmc(4)
      tband=mcmc(5)

      betaw=slice(1)
      gw=slice(2)
      thetaw=slice(3)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0

      do i=1,p
         do j=1,p
            sbinv(i,j)=sb(i,j)
         end do
      end do
      call inverse(sbinv,p,iflag)      

      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()
      
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Slice sampler step for beta, v, and g.
c+++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ Step a  

         call tdbdplogposteri(kk,maxn,nrec,p,
     &                       alpha,beta,theta,gp,v,lambda,
     &                       tau1,tau2,mub,sbinv,x,y,
     &                       workdpw,weight,aux)

         evaly=log(runif())+aux

c+++++++ Step b
         call dbdpstepbslice(0,gp,gw,gL,gR)
         do i=1,(maxn-1)
            do j=1,p
               call dbdpstepbslice(0,beta(i,j),betaw,
     &                             betal(i,j),betar(i,j))
            end do
         end do
         do i=1,maxn
            call dbdpstepbslice(1,theta(i),thetaw,thetal(i),thetar(i))
         end do

c+++++++ Step c
         if(tau1.gt.0.d0)then
            call dbdpstepc1slice(gl,gr,g1)
           else 
            g1=gp
         end if

         do i=1,(maxn-1)
            do j=1,p
               call dbdpstepc1slice(betal(i,j),betar(i,j),beta1(i,j))
            end do
         end do

         do i=1,maxn
            call dbdpstepc1slice(thetal(i),thetar(i),theta1(i))
         end do

         call tdbdplogposteri(kk,maxn,nrec,p,
     &                       alpha,beta1,theta1,g1,v,lambda,
     &                       tau1,tau2,mub,sbinv,x,y,
     &                       workdpw,weight,evalf)

         do while(evaly.gt.evalf)


c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            acrate(1)=acrate(1)+1.d0

            if(tau1.gt.0.d0)then
               call dbdpstepc2slice(0,gp,g1,gl,gr)
               call dbdpstepc1slice(gl,gr,g1)
              else
               g1=gp
            end if
            
            do i=1,(maxn-1)
               do j=1,p
                  call dbdpstepc2slice(0,beta(i,j),beta1(i,j),
     &                                 betal(i,j),betar(i,j))
                  call dbdpstepc1slice(betal(i,j),betar(i,j),
     &                                 beta1(i,j))
               end do
            end do
            do i=1,maxn
               call dbdpstepc2slice(1,theta(i),theta1(i),
     &                              thetal(i),thetar(i))
               call dbdpstepc1slice(thetal(i),thetar(i),theta1(i))
            end do
        
            call tdbdplogposteri(kk,maxn,nrec,p,
     &                          alpha,beta1,theta1,g1,v,lambda,
     &                          tau1,tau2,mub,sbinv,x,y,
     &                          workdpw,weight,evalf)

         end do
        
c+++++++ update parameters

         gp=g1
         do i=1,(maxn-1)
            do j=1,p
               beta(i,j)=beta1(i,j)
            end do
         end do
         do i=1,maxn
            theta(i)=theta1(i)
         end do
 
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ sampling kk
c+++++++++++++++++++++++++++++++++++++++++++++++++
         if(lambda.gt.0.d0)then

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            kk1=kk+1
            if(kk.gt.1)then
               if(runif().gt.0.5)then
                 kk1=kk-1
               end if
            end if

            call tdbdploglike(kk,maxn,nrec,p,beta,theta,x,y,v,
     &                       workdpw,weight,evaly)

            call tdbdploglike(kk1,maxn,nrec,p,beta,theta,x,y,v,
     &                       workdpw,weight,evalf)

            aux=evalf-evaly+
     &          dpoiss(dble(kk1),dble(lambda),1)-
     &          dpoiss(dble(kk) ,dble(lambda),1)

            if(log(dble(runif())).lt.aux)then
               acrate(2)=acrate(2)+1.d0
               kk=kk1
            end if
         end if

c         call intpr("k",-1,kk,1)

c+++++++++++++++++++++++++++++++++++
c+++++++ baseline mean           +++
c+++++++++++++++++++++++++++++++++++

         do i=1,p
            workv1(i)=s0invm(i)
         end do

         do i=1,p
            do j=1,p
               workm1(i,j)=s0inv(i,j)+dble(maxn-1)*sbinv(i,j)
            end do
         end do
         call inverse(workm1,p,iflag)

         do ii=1,(maxn-1)      
            do i=1,p
               workv2(i)=beta(ii,i)
            end do
            workv1=workv1+matmul(sbinv,workv2)
         end do

         workv2=matmul(workm1,workv1)
         call rmvnorm(p,workv2,workm1,workmh1,workv1,mub)

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline covariance matrix +++
c++++++++++++++++++++++++++++++++++++++

         if(tau1.le.0.d0)then
            do i=1,p
               do j=1,p
                  sb(i,j)=psiinv(i,j)
                  workm1(i,j)=0.d0
               end do
               workv1(i)=0.d0
               iflag(i)=0
            end do
  
            do ii=1,(maxn-1)
               do i=1,p
                  do j=1,p
                     sb(i,j)=sb(i,j)+(beta(ii,i)-mub(i))*
     &                               (beta(ii,j)-mub(j))
                  end do
               end do
            end do

            call riwishart(p,nu+maxn-1,sb,sbinv,workm1,workv1,workmh1,
     &                     workmh2,iflag)

         end if

c         call dblepr("sb",-1,sb,p*p)
c         call dblepr("sbinv",-1,sbinv,p*p)

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

               count=0

c+++++++++++++ kk
               count=count+1
               thetasave(isave,count)=kk

c+++++++++++++ mub, sbinv
               do i=1,p
                  count=count+1
                  thetasave(isave,count)=mub(i)
               end do
               do i=1,p
                  do j=i,p
                     count=count+1
                     thetasave(isave,count)=sb(i,j)
                  end do 
               end do

c+++++++++++++ g-prior
               count=count+1
               thetasave(isave,count)=gp

c+++++++++++++ beta
               count=0
               do i=1,(maxn-1)
                  do j=1,p
                     count=count+1
                     randsave(isave,count)=beta(i,j)
                  end do
               end do

c+++++++++++++ stick-brealking weights
               do i=1,maxn
                  count=count+1
                  randsave(isave,count)=theta(i)
               end do

c+++++++++++++ cpo

               do i=1,nrec
                  fw(i)=0.d0
                  do j=1,(maxn-1)
                     aux=0.d0
                     do k=1,p
                        aux=aux+(x(i,k)*beta(j,k))
                     end do
                     v(j)=exp(aux)/(1.d0+exp(aux))
                  end do
                  v(maxn)=1.d0
                  call sickbreak(maxn,v,workdpw,weight)
                  do j=1,maxn
                     call jcomponentbd(theta(j),kk,k)
                     fw(i)=fw(i)+weight(j)*
     &                     dbet(y(i),dble(k),dble(kk-k+1),0)
                  end do
               end do

               do i=1,nrec
                   cpo(i,1)=cpo(i,1)+1.d0/fw(i)
                   cpo(i,2)=cpo(i,2)+fw(i)
               end do 
               
c+++++++++++++ predictions

               do i1=1,npred
                  fm(i1)=0.d0
                  do j=1,(maxn-1)
                     aux=0.d0
                     do k=1,p
                        aux=aux+(xpred(i1,k)*beta(j,k))
                     end do
                     v(j)=exp(aux)/(1.d0+exp(aux))
                  end do
                  v(maxn)=1.d0
                  call sickbreak(maxn,v,workdpw,weight)
                  do j1=1,ngrid
                      fw2(i1,j1)=0.d0
                      do j=1,maxn
                         call jcomponentbd(theta(j),kk,k)
                         fw2(i1,j1)=fw2(i1,j1)+weight(j)*
     &                        dbet(grid(j1),dble(k),dble(kk-k+1),0)
                         if(j1.eq.1)then 
                            fm(i1)=fm(i1)+weight(j)*dble(k)/dble(kk+1)
                         end if
                      end do
                  
                      fmean(i1,j1)=fmean(i1,j1)+fw2(i1,j1)
                  end do

                  meanfpm(i1)=meanfpm(i1)+fm(i1)

                  write(1) (fw2(i1,j1),j1=1,ngrid)
               end do

               write(2) (fm(i),i=1,npred)

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  tmp1=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if
         end if
      end do

      do i1=1,npred
         meanfpm(i1)=meanfpm(i1)/dble(nsave) 
         do j1=1,ngrid
            fmean(i1,j1)=fmean(i1,j1)/dble(nsave)
         end do
      end do 

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)
      end do

      close(unit=1)
      close(unit=2)

      if(cband.eq.1)then
         call hpddensreg(nsave,npred,ngrid,0.05d0,tband,worksam,fs,
     &                   flow,fupp)

         call hpddensregmf(nsave,npred,0.05d0,tband,worksam,
     &                     meanfpl,meanfph)

      end if

      return
      end



c=======================================================================
       subroutine tdbdploglike(kk,maxn,nrec,p,beta,theta,x,y,v,workv,
     &                        weight,eval)
c=======================================================================
c      log-Likelihood for the BDDP 
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,nrec,p
       real(kind=8) beta(maxn-1,p)
       real(kind=8) y(nrec),x(nrec,p),v(maxn)         
       real(kind=8) theta(maxn)

c+++++ External working space
       real(kind=8) workv(maxn+1),weight(maxn)

c+++++ Internal working space
       integer i,j,k
       real(kind=8) dbet,tmp1,tmp2

c+++++ Output
       real(kind=8) eval

c+++++ Algorithm


       eval=0.d0
       do i=1,nrec
          tmp2=0.d0
          do j=1,(maxn-1)
             tmp1=0.d0
             do k=1,p
                tmp1=tmp1+(x(i,k)*beta(j,k))
             end do
             v(j)=exp(tmp1)/(1.d0+exp(tmp1))
          end do
          v(maxn)=1.d0
          call sickbreak(maxn,v,workv,weight)
          do j=1,maxn
             call jcomponentbd(theta(j),kk,k)
             tmp2=tmp2+weight(j)*dbet(y(i),dble(k),dble(kk-k+1),0)
          end do
          eval=eval+log(tmp2)
       end do

       return
       end
 

c=======================================================================
       subroutine tdbdplogpriori(kk,maxn,p,alpha,
     &                           beta,gp,lambda,tau1,tau2,
     &                           theta, mb,sbinv,eval)
c=======================================================================
c      log-priori distribution 
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,p
       real(kind=8) alpha(2),beta(maxn-1,p),gp,lambda,tau1,tau2
       real(kind=8) theta(maxn)
       real(kind=8) mb(p),sbinv(p,p)

c+++++ Output 
       real(kind=8) eval

c+++++ Internal working space
       integer i,j,k
       real(kind=8) dpoiss,dgamm,dbet,tmp1

c+++++ Algorithm

c      kk
       eval=dpoiss(dble(kk),lambda,1)
       eval=eval-log(1.d0-dpoiss(0.d0,lambda,0))

c      g
       if(tau1.gt.0.d0)then
          eval=eval+dgamm(gp,tau1,tau2,1)       
       end if

c      theta
       do i=1,maxn
          eval=eval+dbet(theta(i),alpha(1),alpha(2),1)
       end do

c      beta
       do i=1,(maxn-1)
          tmp1=0.d0
          do j=1,p
             do k=1,p
                tmp1=tmp1+(beta(i,j)-mb(j))*
     &                     sbinv(j,k)*
     &                    (beta(i,k)-mb(k))
             end do
          end do
          eval=eval-0.5d0*tmp1
       end do

       return
       end

c=======================================================================
       subroutine tdbdplogposteri(kk,maxn,nrec,p,
     &                           alpha,beta,theta,gp,v,lambda,
     &                           tau1,tau2,mb,sbinv,x,y,
     &                           workv,weight,eval)
c=======================================================================
c      log-posteriori distribution.
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,nrec,p
       real(kind=8) alpha(2)
       real(kind=8) beta(maxn-1,p),gp,v(maxn),theta(maxn)
       real(kind=8) lambda,tau1,tau2
       real(kind=8) mb(p)
       real(kind=8) sbinv(p,p)
       real(kind=8) x(nrec,p),y(nrec)

c+++++ External working space
       real(kind=8) workv(maxn+1),weight(maxn)

c+++++ Internal working space
       real(kind=8) tmp1

c+++++ Output
       real(kind=8) eval

c+++++ Algorithm

       call tdbdploglike(kk,maxn,nrec,p,beta,theta,x,y,v,workv,
     &                  weight,eval)

       call tdbdplogpriori(kk,maxn,p,
     &                    alpha,beta,gp,lambda,tau1,tau2,theta,                 
     &                    mb,sbinv,tmp1)

       eval=eval+tmp1

       return
       end

