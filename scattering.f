ACHFCSPSE_SCATTERING.  CALCULATING SCATTERING STATES IN LOCAL AND       
1   NON-LOCAL COULOMB-LIKE POTENTIALS.  Z. PAPP.                        
REF. IN COMP. PHYS. COMMUN. 70 (1992) 435                               
*                                                                       
*  program cspse_scattering                                             
*     Zoltan Papp                                                       
*       present address (1991):                                         
*           Inst. Theor. Phys., Univ. Tubingen, Tubingen, Germany       
*       permanent address:                                              
*           Inst. Nucl. Res., Debrecen, Hungary                         
*=================================================================      
*                                                                       
*  potential separable expansion on Coulomb-Sturmian  basis             
*                                                                       
*  Calculating Scattering States in Local and Non-local                 
*  Coulomb-like Potentials                                              
*                                                                       
*=================================================================      
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        character*12 twbp,glag00,glagm05                                
        parameter (ndim=25,npint=48)                                    
        common /iout/ iout,inp                                          
        common /intst/xst(npint),wst(npint),npst                        
        common /intho/xho(npint),who(npint),npho                        
        dimension vp(ndim,ndim)                                         
*                                                                       
        inp=5                                                           
        iout=6                                                          
        read(inp,*) glag00,glagm05                                      
        read(inp,*) xm1,xm2,twbp                                        
        read(inp,*) lm,omcst                                            
        read(inp,*) minn,maxn                                           
        if(maxn.ge.ndim) stop ' maxn > ndim  '                          
        if(minn.gt.maxn) stop ' minn > maxn '                           
*                                                                       
        open(2,file=glag00,status='old')                                
        read(2,*) npst,alfa                                             
        if(alfa.ne.0.0d0) stop ' alfa .ne. 0.0d0 '                      
        read(2,*) (xst(i),wst(i),i=1,npst)                              
        close(2)                                                        
*                                                                       
        open(2,file=glagm05,status='old')                               
        read(2,*) npho,alfa                                             
        if(alfa.ne.-0.5d0) stop ' alfa .ne. -0.5d0 '                    
        read(2,*) (xho(i),who(i),i=1,npho)                              
        close(2)                                                        
*                                                                       
        write(iout,'(/a,a/)') '  Scattering States',                    
     &             ' in Coulomb-like potential'                         
*                                                                       
        xm=xm1*xm2/(xm1+xm2)                                            
        write(iout,'(/1x,a,f12.8,a,f12.8//a,a)') ' mass_1 =',xm1,       
     &    ' mass_2 =',xm2,'  Potential file name : ',twbp               
        write(iout,'(/1x,a,i3,a,f7.4,a,f7.4)') ' lm =',lm,              
     &    '   omcst =',omcst,'   b =',sqrt(xm)*omcst                    
        call cspom(xm,lm,maxn,omcst,twbp,vp,z12,ipflag,1)               
        if(ipflag.lt.0) then                                            
          write(iout,'(a,i2)') ' The potential vanishes for l =',lm     
          stop                                                          
        endif                                                           
        call cstm(xm,lm,z12,omcst,vp,minn,maxn)                         
        stop                                                            
        end                                                             
                                                                        
        subroutine cstm(xm,lm,z12,omcst,vp,minn,maxn)                   
*                                                                       
*                                                                       
*  Formal parameters                                                    
*     xm     - reduced mass                                             
*     lm     - angular momentum                                         
*     omcst  - scalig parameter for the basis                           
*              b=sqrt(xm)*omcst                                         
*     z12    - charge                                                   
*     vp    -  matrix elements of the potential                         
*     minn  -  starting value of the radial quantum number              
*     maxn  -  maximum value of the radial quantum number               
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        complex*16 vv,g,bwf,tl,fl,work,vdet                             
        parameter (ndim=25,hb2=41.801614d0,coul=1.43996518d0,           
     &          pi=3.141592653589793d0)                                 
        common /iout/ iout,inp                                          
        dimension szi(ndim),vp(ndim,ndim),vv(ndim,ndim),g(ndim,ndim),   
     &          bwf(ndim),cst(ndim),kpvt(ndim),work(ndim),vdet(2)       
        szfg(n,i)=(1-exp(-((i-n-1)*6.d0/n)**2))/(1-exp(-6.d0**2))       
*                                                                       
        bst=sqrt(xm)*omcst                                              
*                                                                       
*  <n'l| V_sl|nl>=<n'l| V |nl>-<n'l| V_c |nl>                           
*                                                                       
        do 1 i=1,maxn+1                                                 
           vp(i,i)=vp(i,i)-z12*coul                                     
   1    continue                                                        
*                                                                       
   5    read(inp,*,end=50) en,idr,dr                                    
          write(iout,'(//a,f10.2,/)') '  Initial energy = ',en          
*                                                                       
          do 20 nmm=minn,maxn                                           
            nmax=nmm                                                    
            ind=nmax+1                                                  
*                                                                       
*      calculation of the smoothing coefficients                        
*                                                                       
            do 12 i=1,ind                                               
  12          szi(i)=szfg(ind,i)                                        
            do 13 i=1,ind                                               
              do 13 k=1,i                                               
                 vv(k,i)=szi(i)*szi(k)*vp(k,i)                          
  13        continue                                                    
*                                                                       
*      V^{-1} - G(z)                                                    
*                                                                       
            call zsifa(vv,ndim,ind,kpvt,info)                           
            if(info.ne.0) stop ' singular potential matrix '            
            call zsidi(vv,ndim,ind,kpvt,vdet,work,1)                    
*                                                                       
            call csgrs(en,bst,xm,lm,z12,nmax,g,eta,cst)                 
            do 14 i=1,ind                                               
              do 14 j=1,i                                               
                 vv(j,i)=vv(j,i)-g(j,i)                                 
  14        continue                                                    
*                                                                       
*        cst( ) * ( V^{-1} - G(z) )^{-1} * cst( )                       
*                                                                       
            call zsifa(vv,ndim,ind,kpvt,info)                           
            if(info.ne.0) stop '  singular pse matrix '                 
            do 15 i=1,ind                                               
  15          bwf(i)=cst(i)                                             
            call zsisl(vv,ndim,ind,kpvt,bwf)                            
*                                                                       
            tl=(0.d0,0.d0)                                              
            do 16 i=1,ind                                               
  16          tl=tl+bwf(i)*cst(i)                                       
            zk=sqrt(2*xm/hb2*en)                                        
*                                                                       
* for complex potentials delta should be declaired as complex           
*                                                                       
            delta=log(1-2*(0.d0,1.d0)*zk*tl/en)/(2*(0.d0,1.d0))         
            write(iout,'(/a,i5)') '  nmax=',nmax                        
            write(iout,'(/a,f12.6,a/a,f12.6,a)')                        
     &          '  Coulomb phase shift eta_l   =',eta,' rad ',          
     &          '                      eta_l   =',eta*180/pi,' deg '    
            write(iout,'(/a/8x,a,f12.6,a/8x,a,f12.6,a)')                
     &          '  Colomb modified nuclear ',                           
     &          '  phase shift delta_l =',delta,' rad ',                
     &          '              delta_l =',delta*180/pi,' deg '          
            fl=-tl/en*exp(2*(0.d0,1.d0)*eta)                            
            write(iout,'(/a/a,d15.6,'' +i*'',d15.6)')                   
     &          '  Coulomb modified nuclear',                           
     &          '  scattering amplitude   a''_l = ',fl                  
  20      continue                                                      
*                                                                       
*       calculation of the wave function                                
*                                                                       
          if(idr.ne.0) then                                             
            write(iout,'(/a/)') '  wave function '                      
            call cswfs(en,bst,xm,lm,z12,nmax,bwf,idr,dr,1)              
          endif                                                         
*                                                                       
        goto 5                                                          
  50    return                                                          
        end                                                             
                                                                        
        subroutine cspom(xm,lm,nvmax,omcst,twbp,vp,z12,ipflag,job)      
*                                                                       
*   Calculates the Coulomb-Sturmian matrix elements of the total        
*   V = V_c + V_s  potential.                                           
*                                                                       
*  We suppose that the potential is a fish-bone optical potential       
*  which local part has the form                                        
*                                                                       
*  V_{loc}(r)=v0*(1+v1*r^2)*exp(-beta1*r^2)+v2*exp(-beta2*r^2)          
*      +zch*e^2*erf(beta3*r)/r,                                         
*                                                                       
*  and the nonlocal part is expressed on harmonic oscillator basis.     
*  We suppose also that the potential for high angular momentum         
*  becomes an angular momentum idependent local potential U_{loc}       
*                                                                       
*  The data file 'twbp' is organized as follows.                        
*  iswloc,lkmax  - control parameters: if(iswloc.eq.0) U_{loc}= 0,      
*         maximal angular momentum where V differ from U_{loc}          
*    zchr,v0r,beta1r,v1r,v2r,beta2r,beta3r - parameters for U_{loc}     
*  lk1     - angular momentum                                           
*    zch,v0,beta1,v1,v2,beta2,beta3 - parameters for the local part     
*              of the fish-bone potential                               
*    omfish,npfs - omega parameter of the harmonic oscillator basis     
*      and number of Pauli forbidden states in the fish-bone potential  
*    eta(i) - the first npfs eta(i) stand for the epsilon parameters,   
*             then eta(i) is the eta parameter. if(eta(i).eq.0) then    
*             there are no more non-zero eta parameters.                
*                                                                       
* Formal parameters:                                                    
*   xm     - reduced mass                                               
*   lm     - angular momentum                                           
*   nvmax  - maximal radial quantum number                              
*   omcst  - scalig parameter for the basis                             
*            b=sqrt(xm)*omcst                                           
*   twbp   - data file name for potential parameters                    
*   vp     - potential matrix                                           
*   z12    - charge                                                     
*   ipflag - control parameter                                          
*            ipflag=-2    V vanishes for lm                             
*            ipflag=-1    V vanishes for lm                             
*            ipflag= 0    V angular momentum independent                
*            ipflag= 1    V depends on the angular momentum lm          
*   job    - control parameter                                          
*            job=0    V-U_{loc} is calculated                           
*            job=1    V is calculated                                   
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        real*8 mbar                                                     
        character*12 twbp                                               
        parameter (ndim=25,npint=48)                                    
        parameter (hb2=41.801614d0)                                     
*        parameter (hb2=1.0d0)                                          
        common /potpar/ zch,v0,v1,beta1,v2,beta2,beta3                  
        common /potpr/ zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r,iswloc     
        common /intho/ xho(npint),who(npint),npho                       
        common /intst/ xst(npint),wst(npint),npst                       
        dimension ho(ndim,npint),st(ndim,npint),eta(ndim),              
     &    mbar(ndim,ndim),tkii(ndim),tkmi(ndim),tr(ndim,ndim),          
     &    vfb(ndim,ndim),vp(ndim,ndim),pf(npint)                        
*                                                                       
* reading in the potential parameters                                   
*                                                                       
        open(2,file=twbp,status='old')                                  
        read(2, *) iswloc,lkmax                                         
        read(2, *) zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r                
*                                                                       
   1    read(2, *,end=4) lk1                                            
        read(2, *) zch,v0,v1,beta1,v2,beta2,beta3                       
*                                                                       
        read(2,*) omfish,npfs                                           
        i=0                                                             
   2    i=i+1                                                           
        read(2,*) eta(i)                                                
        if(i.le.npfs.or.eta(i).ne.0.d0) goto 2                          
        neta=i-1                                                        
        if(lk1.ne.lm) then                                              
          goto 1                                                        
        else                                                            
          close(2)                                                      
          ipflag=2                                                      
          z12=zch                                                       
          goto 5                                                        
        endif                                                           
   4    close(2)                                                        
        if(job.eq.1) then                                               
          if(iswloc.ne.0) then                                          
            ipflag=1                                                    
            z12=zchr                                                    
            neta=0                                                      
          else                                                          
            ipflag=0                                                    
            return                                                      
          endif                                                         
        elseif(job.eq.0) then                                           
          if(lm.lt.lkmax) then                                          
            ipflag=-1                                                   
          else                                                          
            ipflag=-2                                                   
          endif                                                         
          return                                                        
        else                                                            
          stop  ' error in calling cspom '                              
        endif                                                           
   5    continue                                                        
*                                                                       
*  calculation of the matrix elements of the local part of the          
*  potential on Coulomb-Sturmian basis                                  
*                                                                       
        nn=nvmax+1                                                      
        bst=sqrt(xm)*omcst                                              
        do 6 i=1,npst                                                   
          call csturm(st(1,i),xst(i),lm,nvmax)                          
          xib=xst(i)/(2*bst)                                            
          if(job.eq.1) then                                             
            if(ipflag.eq.1) then                                        
              pf(i)=potfr(xib)/(2*bst)                                  
            else                                                        
              pf(i)=potfun(xib)/(2*bst)                                 
            endif                                                       
          else                                                          
            if(iswloc.eq.0) then                                        
              pf(i)=potfun(xib)/(2*bst)                                 
            else                                                        
              pf(i)=(potfun(xib)-potfr(xib))/(2*bst)                    
            endif                                                       
          endif                                                         
   6    continue                                                        
        do 8 i=1,nn                                                     
          do 8 j=1,i                                                    
            vv=0.0d0                                                    
            do 7 n=1,npst                                               
   7           vv=vv+st(i,n)*st(j,n)*pf(n)*wst(n)                       
            vp(j,i)=vv                                                  
            vp(i,j)=vv                                                  
   8    continue                                                        
*                                                                       
        if(neta.ne.0) then                                              
*                                                                       
*  calculation of the matrix elements of the non-local part of the      
*  potential on Coulomb-Sturmian basis                                  
*                                                                       
*   calculate the norm kernel eigenfunctions                            
*   (harmonic oscillator functions in tis case)                         
*                                                                       
          bho=1.d0/sqrt(xm*omfish)                                      
          do 9 i=1,npho                                                 
            xsq=sqrt(xho(i))                                            
            call hof(ho(1,i),xsq,lm,nvmax)                              
            pf(i)=potfun(xsq*bho)                                       
   9      continue                                                      
          do 10 i=neta+1,nn                                             
  10        eta(i)=0.d0                                                 
*                                                                       
*   calculate the M-bar matrix and the matrix elements of the           
*   kinetic energy                                                      
*                                                                       
          hbm=hb2/2*omfish                                              
          do 12 j=1,nn                                                  
          do 12 i=1,j                                                   
          if(i.le.npfs) then                                            
             mbar(i,j)=1.d0                                             
          else                                                          
            if(i.eq.j.or.i.gt.neta) then                                
             mbar(i,j)=0.d0                                             
            else                                                        
             mbar(i,j)=1.d0-sqrt((1.d0-eta(i))/(1.d0-eta(j)))           
            endif                                                       
          endif                                                         
          if(i.le.neta) then                                            
            if(i.eq.j) then                                             
              tkii(i)=hbm*(2*i+lm-0.5d0)                                
              if(i.le.npfs) tkii(i)=tkii(i)-eta(i)                      
            elseif(i.eq.j-1) then                                       
              tkmi(i)=hbm*sqrt(i*(i+lm+0.5d0))                          
            endif                                                       
          endif                                                         
  12      continue                                                      
*                                                                       
*   calcullate the harmonic oscillator matrix elements of the           
*   non-local part of the potential                                     
*                                                                       
          do 17 j=1,nn                                                  
          do 17 i=1,j                                                   
          if(mbar(i,j).eq.0.d0) then                                    
            vfb(i,j)=0.0d0                                              
          else                                                          
            vv=0.0d0                                                    
            do 15 n=1,npho                                              
  15        vv=vv+ho(i,n)*ho(j,n)*pf(n)*who(n)                          
            if(i.eq.j) then                                             
              vv=vv+tkii(i)                                             
            elseif(i.eq.j-1) then                                       
              vv=vv+tkmi(i)                                             
            endif                                                       
            vfb(i,j)=vv*mbar(i,j)                                       
          endif                                                         
          vfb(j,i)=vfb(i,j)                                             
  17      continue                                                      
*                                                                       
*  transform the matrix elemnts from harmonic oscillator basis          
*  to Coulomb-Sturmian basis                                            
*                                                                       
          call csho(lm,nvmax,bst,bho,tr)                                
*                                                                       
          do 20 j=1,nn                                                  
          do 20 i=1,j                                                   
            vv=0.0d0                                                    
            do 19 k1=1,nn                                               
            do 19 k2=1,nn                                               
               vv=vv+tr(i,k1)*vfb(k1,k2)*tr(j,k2)                       
  19        continue                                                    
            vp(i,j)=vp(i,j)-vv                                          
            vp(j,i)=vp(i,j)                                             
  20      continue                                                      
        endif                                                           
*                                                                       
        return                                                          
        end                                                             
                                                                        
        subroutine csho(l,nvmax,bst,bho,tr)                             
*                                                                       
*   Calculates the overlap matrix elements betwen Coulomb-Sturmian      
*   and norm kernel eigenfunctions (harmonic  oscillator functions)     
*                                                                       
*   l     - angular momentum                                            
*   nvmax - maximal radial quantum number                               
*   bst   - b parameter of the Coulomb-Sturmian functions               
*   bho   - b parameter of the harmonic oscillator functions            
*   tr    - overlap metrix elements                                     
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        parameter (ndim=25,npint=48)                                    
        common /intst/xst(npint),wst(npint),npst                        
        dimension  h(ndim,npint),s(ndim,npint),tr(ndim,ndim)            
        do 1 i=1,npst                                                   
        call csturm(s(1,i),2*xst(i),l,nvmax)                            
        xib=xst(i)/bst                                                  
        call howf(h(1,i),xib,bho,l,nvmax)                               
   1    continue                                                        
        do 8 ni=1,nvmax+1                                               
        do 8 nj=1,nvmax+1                                               
        vv=0.d0                                                         
        do 3 i=1,npst                                                   
   3    vv=vv+wst(i)*s(ni,i)*h(nj,i)                                    
        tr(ni,nj)=vv/bst                                                
   8    continue                                                        
        return                                                          
        end                                                             
                                                                        
                                                                        
        subroutine csturm(st,x,l,nmax)                                  
*                                                                       
*   Calculates Coulomb-Sturmian fuctions                                
*   (without the exponential factor)                                    
*                                                                       
*   st(n+1)=[n!/(n+2l+1)!]^{1/2} x^{l+1} L_n^{2l+1} (x)                 
*           from n=0 to nmax                                            
*   x  - radial variable (2b*r)                                         
*   l  - angular momentum                                               
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        dimension st(1)                                                 
        l2=2*l+1                                                        
        sfl=1.d0                                                        
        do 1 i=2,l2                                                     
   1    sfl=sfl*i                                                       
        xl0=0.d0                                                        
        st(1)=x**(l+1)/sqrt(sfl)                                        
        do 5 i=1,nmax                                                   
        df=sqrt(dble(i*(i+l2)))                                         
        st(i+1)=((2*i+l2-1-x)*st(i)-xl0)/df                             
        xl0=df*st(i)                                                    
   5    continue                                                        
        return                                                          
        end                                                             
                                                                        
        subroutine hof(ho,x,l,nmax)                                     
*                                                                       
*  ho(n+1)=howf(n,l,r,r0)*sqrt(r0/2)*exp((r/r0)**2/2)                   
*  where howf is the harmonic oscillator function with                  
*  parameters                                                           
*   n:  radial quantum number                                           
*   l:  angular momentum variable                                       
*   r:  radial variable                                                 
*   r0: oscillator lenght parameter                                     
*  --Other input parameters--                                           
*   x:  r/r0                                                            
*   nn: maximal radial quantum number                                   
*  --Output parameters--                                                
*   ho: howf(n,l,r,r0)*sqrt(r0/2)*exp((r/r0)**2/2)                      
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        parameter (sqpi=1.772453850905516d0)                            
        dimension ho(1)                                                 
        x2=x*x                                                          
        cl=2*x2/sqpi                                                    
        do 1 i=1,l                                                      
    1   cl=2*cl/(2*i+1)*x2                                              
        ho(1)=sqrt(cl)                                                  
        rn1=0.d0                                                        
        do 2 i=1,nmax                                                   
        df=sqrt(i*(i+l+0.5d0))                                          
        ho(i+1)=((2*i+l-0.5d0-x2)*ho(i)-rn1)/df                         
        rn1=df*ho(i)                                                    
    2   continue                                                        
        return                                                          
        end                                                             
                                                                        
        subroutine howf(ho,r,r0,l,nmax)                                 
*                                                                       
*  ho(n+1)=howf(n,l,r,r0)                                               
*  where howf is the harmonic oscillator function with                  
*  parameters                                                           
*    n:  radial quantum number                                          
*    l:  angular momentum variable                                      
*    r:  radial variable                                                
*    r0: oscillator lenght parameter                                    
*  --Other input parameters--                                           
*  nn: maximal radial quantum number                                    
*  --Output parameters--                                                
*  ho: howf(n,l,r,r0)                                                   
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        dimension ho(1)                                                 
        x=r/r0                                                          
        aa=sqrt(2/r0)*exp(-x**2/2)                                      
        call hof(ho,x,l,nmax)                                           
        do 1 i=0,nmax                                                   
   1    ho(i+1)=aa*ho(i+1)                                              
        return                                                          
        end                                                             
                                                                        
        function potfun(r)                                              
*                                                                       
*  The potential parameters are stoored in the commonn block /potpar/.  
*  In this case they describe the following potential:                  
*  rpot=v0*(1+v1*r^2)*exp(-beta1*r^2)+v2*exp(-beta2*r^2)                
*      +zch*coul*erf(beta3*r)/r                                         
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        common /potpar/ zch,v0,v1,beta1,v2,beta2,beta3                  
        common /potpr/ zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r,iswloc     
        parameter (coul=1.43996518d0)                                   
*        parameter (coul=1.0d0)                                         
*                                                                       
        rr=r*r                                                          
        rpot=v0*(1.d0+v1*rr)*exp(-beta1*rr)                             
     &          +v2*exp(-beta2*rr)                                      
        if(zch.eq.0.d0) then                                            
          rpotc=0.d0                                                    
        else                                                            
          rpotc=zch*coul*erf(beta3*r)/r                                 
        endif                                                           
        potfun=rpot+rpotc                                               
        return                                                          
*                                                                       
        entry potfr(r)                                                  
*                                                                       
*       local auxiliary potential U_{loc}                               
*                                                                       
*  The potential parameters are stoored in the commonn block /potpr/.   
*  In this case they describe the following potential:                  
*  rpot=v0r*(1+v1r*r^2)*exp(-beta1r*r^2)+v2r*exp(-beta2r*r^2)           
*      +zchr*coul*erf(beta3r*r)/r                                       
*                                                                       
        rr=r*r                                                          
        rpot=v0r*(1.d0+v1r*rr)*exp(-beta1r*rr)                          
     &          +v2r*exp(-beta2r*rr)                                    
        if(zchr.eq.0.d0) then                                           
          rpotc=0.d0                                                    
        else                                                            
          rpotc=zchr*coul*erf(beta3r*r)/r                               
        endif                                                           
        potfr=rpot+rpotc                                                
        return                                                          
        end                                                             
                                                                        
        function erf(x)                                                 
*                                                                       
*       calculates the error function  (accuracy > 1.2d-7 )             
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        z=abs(x)                                                        
        t=1.d0/(1.d0+0.5d0*z)                                           
        erfc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+                   
     &    t*(0.37409196d0+t*(0.09678418d0+t*(-0.18628806d0+             
     &    t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0+             
     &    t*(-0.82215223d0+t*0.17087277)))))))))                        
        if(x.lt.0.d0) erfc=2.d0-erfc                                    
        erf=1.d0-erfc                                                   
        return                                                          
        end                                                             
                                                                        
        subroutine csgrs(e,bst,xm,l,z12,nmax,g,eta,cst)                 
*                                                                       
*       Calculates the matrix elements of the                           
*       Coulomb-Green operator on Coulomb-Sturmian basis                
*                                                                       
*     e      - energy parameter of the Green operator                   
*     bst    - b parameter for the basis                                
*     xm     - reduced mass                                             
*     lm     - angular momentum                                         
*     z12    - charge                                                   
*     nmax   - maximum value of the radial quantum number               
*     g      - matrix elements of the Green operator                    
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        complex*16 g,zki,gammai,s0,f21,clogam                           
        parameter (ndim=25,hb2=41.801614d0,coul=1.43996518d0,           
     &          pi=3.141592653589793d0)                                 
        dimension g(ndim,ndim),cst(ndim)                                
*                                                                       
        zk=sqrt(2*xm/hb2*e)                                             
        zki=(0.0d0,1.0d0)*zk                                            
        ctr=z12*coul*xm/hb2                                             
        gamma=ctr/zk                                                    
        gammai=(0.0d0,1.0d0)*gamma                                      
*                                                                       
*  < 0l | \Delta g_l^c (E) \Delta | 0l>                                 
*                                                                       
        g(1,1)=-4*xm*bst/(hb2*(l+gammai+1)*(bst-zki)**2)                
     &          *f21(-l+gammai,l+gammai+2,((bst+zki)/(bst-zki))**2)     
*                                                                       
        zkb2=zk**2+bst**2                                               
        sx1=(zk**2-bst**2)/zkb2                                         
        sx2=4*ctr*bst/zkb2                                              
        s0=4*xm/hb2*bst/zkb2                                            
*                                                                       
*  < 0l | \Delta g_l^c (E) \Delta | il>                                 
*                                                                       
        do 5 i=0,nmax-1                                                 
          df=sqrt(dble((i+1)*(i+2*l+2)))                                
          g(1,i+2)=(g(1,i+1)*(2*(i+l+1)*sx1-sx2)-s0)/df                 
          g(i+2,1)=g(1,i+2)                                             
          s0=g(1,i+1)*df                                                
   5    continue                                                        
*                                                                       
*  < jl | \Delta g_l^c (E) \Delta | il>                                 
*                                                                       
        do 6 i=0,nmax-1                                                 
          s0=(0.d0,0.d0)                                                
          do 6 j=0,i                                                    
            df=sqrt(dble((j+1)*(j+2*l+2)))                              
            g(j+2,i+2)=(g(j+1,i+2)*(2*(j+l+1)*sx1-sx2)-s0)/df           
            g(i+2,j+2)=g(j+2,i+2)                                       
            s0=g(j+1,i+2)*df                                            
   6    continue                                                        
*                                                                       
        zb=zk/bst                                                       
        if(z12.eq.0.0d0) then                                           
          eta=0.d0                                                      
          cst(1)=(2*zb/(1+zb*zb))**(l+1)                                
        else                                                            
          s0=clogam(l+gammai+1)                                         
          eta=dimag(s0)                                                 
          cst(1)=sqrt(2*pi*gamma/(exp(2*pi*gamma)-1))*                  
     &           exp(2*gamma*atan(zb))*(2*zb/(1+zb*zb))**(l+1)          
        endif                                                           
        do 7 i=1,l                                                      
   7      cst(1)=cst(1)*sqrt((gamma*gamma+i*i)/(i*(i+0.5d0)))           
        s2=0.d0                                                         
        do 8 i=0,nmax-1                                                 
          df=sqrt(dble((i+1)*(i+2*l+2)))                                
          cst(i+2)=(cst(i+1)*(2*(i+l+1)*sx1-sx2)-s2)/df                 
          s2=cst(i+1)*df                                                
   8    continue                                                        
*                                                                       
        return                                                          
        end                                                             
                                                                        
        subroutine cswfs(e,bst,xm,l,z12,nmax,bwf,idr,dr,isw)            
*                                                                       
*       Calculates the wave function                                    
*                                                                       
*     e      - energy parameter of the Coulomb-Green operator           
*     bst    - b parameter for the basis                                
*     xm     - reduced mass                                             
*     lm     - angular momentum                                         
*     z12    - charge                                                   
*     nmax   - maximum value of the radial quantum number               
*     bwf    - wave function coefficients b                             
*     idr    - number of points in calculating wave function            
*     dr     - step size in calculating wave function                   
*     isw    - if(isw.ne.0) the Coulomb function is added on            
*                                                                       
        implicit real*8 (a-h,o-z)                                       
        complex*16 bwf,zki,gammai,xum,cheby,gcs,gexpi,a0,s              
        parameter (ndim=25,ldim=100,eps=1.d-15,step=100)                
        parameter (hb2=41.801614d0,coul=1.43996518d0)                   
        common /iout/ iout                                              
        dimension bwf(ndim),cheby(ldim),gcs(ndim),                      
     &    fc(2*ndim),fcp(1),gc(1),gcp(1)                                
        equivalence (gcs(1),fc(1))                                      
*                                                                       
        sfl=1                                                           
        do 1 i=1,2*l+1                                                  
   1      sfl=sfl*i                                                     
        sfl=(2*bst)**(l+1)/sqrt(sfl)                                    
        zk=sqrt(2*xm/hb2*e)                                             
        zki=(0.d0,1.d0)*zk                                              
        ctr=z12*coul*xm/hb2                                             
        gamma=ctr/zk                                                    
        gammai=(0.0d0,1.0d0)*gamma                                      
*                                                                       
        zkb2=zk**2+bst**2                                               
        sx1=(zk**2-bst**2)/zkb2                                         
        sx2=4*ctr*bst/zkb2                                              
        s0=4*xm/hb2*bst/zkb2                                            
*                                                                       
        xum=-2*xm/hb2*sfl/(bst-zki)*                                    
     &          exp((l-gammai)*log(2*zki/(zki-bst)))                    
        call ch1f0(-l+gammai,-(bst+zki)/(2*zki),ldim,cheby)             
        nchu=10                                                         
        do 8 np=1,idr                                                   
          r=np*dr                                                       
          rho=zk*r                                                      
*                                                                       
*  < r |  g_l^c (E) \Delta | 0l>                                        
*                                                                       
          gcs(1)=xum*r**(l+1)*exp(-bst*r)*                              
     &          gexpi(-(bst+zki)*r,l+gammai,cheby,nchu)                 
          st0=0.d0                                                      
          df=0.d0                                                       
          st1=sfl*exp(-bst*r)*r**(l+1)                                  
          a0=(0.d0,0.d0)                                                
*                                                                       
*  < r |  g_l^c (E) \Delta | il>                                        
*                                                                       
          do 4 i=0,nmax-1                                               
            st2=2*(i+l+1-bst*r)*st1-df*st0                              
            df=sqrt(dble((i+1)*(i+2*l+2)))                              
            gcs(i+2)=(gcs(i+1)*(2*(i+l+1)*sx1-sx2)-a0-s0*st1)/df        
            a0=df*gcs(i+1)                                              
            st0=st1                                                     
            st1=st2/df                                                  
   4      continue                                                      
*                                                                       
*  < r| u_l>=\sum_{i=0}^{nmax} bwf(i+1) * <r|g_l^c (E) \Delta|il>       
*                                                                       
          s=(0.d0,0.d0)                                                 
          do 5 i=1,nmax+1                                               
   5        s=s+bwf(i)*gcs(i)                                           
*                                                                       
          if(isw.ne.0) then                                             
             call rcwff(rho,gamma,l,l,fc,gc,fcp,gcp,eps,step,1)         
             s=fc(l+1)+s                                                
          endif                                                         
*                                                                       
          write(iout,'(2x,f6.2,2d12.4)')  r,s                           
*                                                                       
   8    continue                                                        
*                                                                       
        return                                                          
        end                                                             
                                                                        
        function f21(a,c,z)                                             
*                                                                       
*  calculates the hypergeometric function 2_F_1 ( a, 1; c; z )          
*                                                                       
        complex*16 f21,a,c,z,p0,q0,q1,r,p2,q2,qq                        
        real*8 eps                                                      
        parameter (eps=1.d-14,nn=200)                                   
        common /iout/ iout                                              
*                                                                       
        p0=(1.d0,0.d0)                                                  
        q0=(0.d0,0.d0)                                                  
        q1=(1.d0,0.d0)                                                  
        qq=(0.d0,0.d0)                                                  
        n=0                                                             
   1    if(abs(1-qq/q1).gt.eps) then                                    
           r=-(a+n)*(c+n-1)/((c+2*n-1)*(c+2*n))*z                       
           p2=1+r*p0                                                    
           q2=q1+r*q0                                                   
           p0=1/p2                                                      
           q0=q1*p0                                                     
           qq=q2*p0                                                     
           n=n+1                                                        
           r=-n*(c+n-a-1)/((c+2*n-2)*(c+2*n-1))*z                       
           p2=1+r*p0                                                    
           q2=qq+r*q0                                                   
           p0=1/p2                                                      
           q0=qq*p0                                                     
           q1=q2*p0                                                     
           if(n.gt.nn) then                                             
              f21=q1                                                    
              write(iout,'('' estimated accuracy in f21 '',d16.9)')     
     &          abs((qq-q1)/q1)                                         
              return                                                    
           endif                                                        
           goto 1                                                       
        else                                                            
           f21=q1                                                       
           return                                                       
        endif                                                           
        end                                                             
                                                                        
        function clogam(z)                                              
c                                                                       
c    computes the logarithm of the gamma function                       
c    for any complex argument K.S Kolbig, CPC 4 (1972) 221              
c                                                                       
        implicit real*8 (a-h,o-z)                                       
        complex*16 z,v,h,r,clogam                                       
        common /iout/ iout                                              
        dimension b(10)                                                 
c                                                                       
        data lerr /3/                                                   
        data pi /3.14159 26535 89793d0/                                 
        data b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10)         
     &       /+8.33333 33333 3333d-2, -2.77777 77777 7778d-3,           
     &        +7.93650 79365 0794d-4, -5.95238 09523 8095d-4,           
     &        +8.41750 84175 0842d-4, -1.91752 69175 2692d-3,           
     &        +6.41025 64102 5641d-3, -2.95506 53594 7712d-2,           
     &        +1.79644 37236 8831d-1, -1.39243 22169 0590d+0/           
c                                                                       
        x=dble(z)                                                       
        t=dimag(z)                                                      
        if(-abs(x) .eq. aint(x) .and. t .eq. 0.d0) go to 5              
        f=abs(t)                                                        
        v=dcmplx(x,f)                                                   
        if(x .lt. 0.d0) v=1.d0-v                                        
        h=0.d0                                                          
        c=dble(v)                                                       
        if(c .ge. 7.d0) go to 3                                         
        n=6-int(c)                                                      
        h=v                                                             
        d=dimag(v)                                                      
        a=atan2(d,c)                                                    
        if(n .eq. 0) go to 2                                            
        do 1 i = 1,n                                                    
        c=c+1.d0                                                        
        v=dcmplx(c,d)                                                   
        h=h*v                                                           
   1    a=a+atan2(d,c)                                                  
   2    h=dcmplx(0.5d0*log(dble(h)**2+dimag(h)**2),a)                   
        v=v+1.d0                                                        
   3    r=1.d0/v**2                                                     
        clogam=.918938533204673d0+(v-.5d0)*log(v)-v+(b(1)+r*(b(2)+      
     &  r*(b(3)+r*(b(4)+r*(b(5)+r*(b(6)+r*(b(7)+r*(b(8)+r*(b(9)+        
     &  r*b(10))))))))))/v-h                                            
        if(x .ge. 0.d0) go to 4                                         
c                                                                       
        a=aint(x)-1.d0                                                  
        c=pi*(x-a)                                                      
        d=pi*f                                                          
        e=exp(-2.d0*d)                                                  
        f=sin(c)                                                        
        e=d+0.5d0*log(e*f**2+0.25d0*(1.d0-e)**2)                        
        f=atan2(cos(c)*tanh(d),f)-a*pi                                  
        clogam=1.1447298858494d0-dcmplx(e,f)-clogam                     
c                                                                       
   4    if(sign(1.d0,t) .eq. -1.d0) clogam=conjg(clogam)                
        return                                                          
c                                                                       
   5    write(iout,100) x                                               
        clogam=(0.d0,0.d0)                                              
        return                                                          
 100    format(' clogam ... argument is non positive integer = ',f20.2) 
c                                                                       
        end                                                             
                                                                        
        function gexpi(a,u,b,n)                                         
*                                                                       
*  calculates                                                           
*  gexpi=exp(-a) * integral from 0 to 1 {exp(a*t) * t**u * f(t)} dt ,   
*  where f(t) is a binomial function, expanded in terms of              
*  Chebyshev polinomials, f(t) = sum_{i=0}^n b_i T^{*}_i.               
*  See Y.L. Luke: Algorithms for the Computation of Mathematical        
*  Functions, Academic Press, New York, 1977, p. 126-129,               
*  formula (19) and (20).                                               
*                                                                       
        implicit complex*16 (a-h,o-z)                                   
        real*8 eps,test                                                 
        parameter (ldim=100,ml=1,mu=2,m=4,nband=5,eps=1.d-4)            
        common /iout/ iout                                              
        dimension b(ldim),g(ldim),band(nband,ldim),ipvt(ldim)           
*                                                                       
*  set up the system of linear band equations                           
*  for the Chebyshev coefficients of gexpi                              
*                                                                       
        n=min(n,ldim)                                                   
*                                                                       
   1    g(1)=2*b(1)-b(2)                                                
        do 2 i=2,n-1                                                    
          g(i)=b(i)-b(i+1)                                              
   2    continue                                                        
        g(n)=b(n)                                                       
        do 5 j=1,n                                                      
          i1=max(1,j-mu)                                                
          i2=min(n,j+ml)                                                
          do 3 i=i1,i2                                                  
            k=i-j+m                                                     
            if(i.eq.j+1) then                                           
              band(k,j)=a/4                                             
            elseif(i.eq.j) then                                         
              band(k,j)=a/4+u+i                                         
            elseif(i.eq.j-1) then                                       
              band(k,j)=-a/4-u+i-1                                      
              if(i.eq.1) band(k,j)=-u                                   
            elseif(i.eq.j-2) then                                       
              band(k,j)=-a/4                                            
            endif                                                       
   3      continue                                                      
   5    continue                                                        
*                                                                       
*  solution of the system linear equations                              
*                                                                       
        call zgbfa(band,nband,n,ml,mu,ipvt,info)                        
        if(info.ne.0) write(iout,'('' zgbfa info ='',i3)') info         
        call zgbsl(band,nband,n,ml,mu,ipvt,g,0)                         
        g(1)=g(1)/2                                                     
*                                                                       
*  gexpi=sum_{i=0}^{n-1} g(i), 1/(u+1)=sum_{i=0}^{n-1} (-1)**i*g(i)     
*                                                                       
        f0=0.d0                                                         
        ip=1                                                            
        ss=(0.d0,0.d0)                                                  
        do 8 i=1,n                                                      
          f0=f0+ip*g(i)                                                 
          ip=-ip                                                        
          ss=ss+g(i)                                                    
   8    continue                                                        
        test=abs(f0-1/(u+1))                                            
        if(test.gt.eps) then                                            
          if(n.eq.ldim) then                                            
            write(iout,*) ' estimated accuracy in gexpi ', test         
            goto 9                                                      
          endif                                                         
          n=min(n+10,ldim)                                              
          goto 1                                                        
        endif                                                           
   9    gexpi=ss                                                        
*                                                                       
        return                                                          
        end                                                             
                                                                        
        subroutine ch1f0(ap,w,nm,c)                                     
*                                                                       
*  Calculates the Chebyshev expansion of the binomial function.         
*  See Y.L. Luke: Algorithms for the Computation of Mathematical        
*  Functions, Academic Press, New York, 1977, p. 53.                    
*                                                                       
*       (1+w*t)^(-ap)=sum_{i=0}^{nm-1} c_i (w) * T^{*}_i (t)            
*                                                                       
        implicit complex*16 (a-h,o-z)                                   
        real*8 eps,test1                                                
        parameter (eps=1.d-10)                                          
        common /iout/ iout                                              
        dimension c(1)                                                  
*                                                                       
*  test wheather ap is a negative integer                               
*                                                                       
        iap=iabs(int(dble(ap)))                                         
        if(ap.eq.dcmplx(-iap)) then                                     
          n=iap-1                                                       
          do 1 i=n+2,nm                                                 
   1         c(i)=0.d0                                                  
        else                                                            
          n=nm-2                                                        
        endif                                                           
*                                                                       
*  computes coefficients by means of backward recurence scheme          
*                                                                       
        start=(2.d0*(1.d0+2.d0/dble(w)))**(-n)                          
        a2=0.d0                                                         
        a1=start                                                        
        nc=n+2                                                          
        c(nc)=start                                                     
        x1=n+2                                                          
        v1=1.d0-ap                                                      
        afac=2.d0+4.d0/w                                                
        do 2 k=1,n+1                                                    
           x1=x1-1.d0                                                   
           nc=nc-1                                                      
           c(nc)=-(x1*afac*a1+(x1+v1)*a2)/(x1-v1)                       
           a2=a1                                                        
           a1=c(nc)                                                     
   2    continue                                                        
        c(1)=c(1)/2.d0                                                  
        rho=c(1)                                                        
        p=1.d0                                                          
        do 3 i=2,n+2                                                    
           rho=rho-p*c(i)                                               
           p=-p                                                         
   3    continue                                                        
        do 4 i=1,n+2                                                    
           c(i)=c(i)/rho                                                
   4    continue                                                        
*                                                                       
*   accuracy test                                                       
*                                                                       
        if(n.gt.1) then                                                 
           sum1=-c(2)+4.d0*c(3)                                         
           p=-1.d0                                                      
           do 5 i=2,n                                                   
              sum1=sum1+p*(i+1)*(i+1)*c(i+2)                            
              p=-p                                                      
   5       continue                                                     
           test1=abs(sum1-ap*w/2.d0)                                    
           if(test1.gt.eps)  write(iout,'(a,d12.3)')                    
     &        ' estimated accuracy in ch1f0 : ',test1                   
        endif                                                           
*                                                                       
        return                                                          
        end                                                             
*                                                                       
*  further subroutines from the LINPACK library:                        
*     zsifa, zsidi, zsisl, zgbfa, zgbsl, zgbcal, dcabs1,                
*     izamax, zdotu, zaxpy, zdotc    and                                
*  Coulomb function RCWF by                                             
*   A.R.Barnett, D.H.Feng, J.W.Steed,L.J.B.Goldfarb, CPC 8 (1974) 377   
*--------------------------------------------------------------------   
*                                                                       
*  one may run the code in UNIX as  a.out < ao.sd > ao.sr               
*                                                                       
*  input file  ao.sd                                                    
*                                                                       
'gisp00.d','gispm05.d',                                                 
 4.0235d0,16.088d0,'ao16p.scd',                                         
 3, 1.75,                                                               
 24, 24,                                                                
 12.d0, 25, .2,                                                         
*                                                                       
*  input file gisp00.d                                                  
*                                                                       
          48                 0.                                         
  2.981123582996164E-002  7.426200582802977E-002                        
  0.157107990617875       0.152271949809351                             
  0.386265037576454       0.190409088263911                             
  0.717574694116971       0.186633059484807                             
   1.15139383402644       0.153424200157578                             
   1.68818582341905       0.108779692807490                             
   2.32852700665323       6.746073860921927E-002                        
   3.07311086165265       3.688119411582105E-002                        
   3.92275241304649       1.785684426915667E-002                        
   4.87839335592135       7.677616514497570E-003                        
   5.94110805462456       2.935785903739469E-003                        
   7.11211053589075       9.990655378158767E-004                        
   8.39276259909124       3.025980169922555E-004                        
   9.78458318468734       8.153871180355333E-005                        
   11.2892591680095       1.953158715728041E-005                        
   12.9086577782855       4.154182945052104E-006                        
   14.6448408832097       7.833700380277454E-007                        
   16.5000814289646       1.307394774920574E-007                        
   18.4768823868741       1.927071408017016E-008                        
   20.5779986340222       2.502638937126329E-009                        
   22.8064622905214       2.855785508771569E-010                        
   25.1656121564391       2.854622412059074E-011                        
   27.6591280444806       2.491010684937161E-012                        
   30.2910710010086       1.890336606971472E-013                        
   33.0659306624988       1.242162685949120E-014                        
   35.9886813274790       7.034231520212479E-016                        
   39.0648487641978       3.414549148591794E-017                        
   42.3005903629031       1.412315414895672E-018                        
   45.7027920385115       4.944218008097224E-020                        
   49.2791863828369       1.453952481367813E-021                        
   53.0384980878167       3.561068365003849E-023                        
   56.9906248148045       7.194055996494519E-025                        
   61.1468647861403       1.185537228350548E-026                        
   65.5202069290186       1.573491357075550E-028                        
   70.1257062361133       1.657285440919349E-030                        
   74.9809775189113       1.361434162716336E-032                        
   80.1068573503244       8.546155813962816E-035                        
   85.5283111160343       4.000090532480939E-037                        
   91.2757079936682       1.355019991102786E-039                        
   97.3866677135817       3.201636795354339E-042                        
   103.908833357176       5.035869166060774E-045                        
   110.904220884976       4.962487540702082E-048                        
   118.456425046284       2.823510716119612E-051                        
   126.683425768886       8.268446069503012E-055                        
   135.762589577865       1.049064847820770E-058                        
   145.986432709463       4.346574422738855E-063                        
   157.915612022978       3.434736438395980E-068                        
   172.996328148563       1.319066088398003E-074                        
*                                                                       
*  input file gispm05.d                                                 
*                                                                       
          48 -0.500000000000000                                         
  1.278457233711745E-002  0.446540046995086                             
  0.115081483588378       0.403226026958469                             
  0.319783878279432       0.328759261196996                             
  0.627109535856757       0.241966233725031                             
   1.03738672215406       0.160707900161284                             
   1.55105613079625       9.627981346217995E-002                        
   2.16867351479483       5.200068054248895E-002                        
   2.89091304669443       2.530275022560015E-002                        
   3.71857145710970       1.108326257887788E-002                        
   4.65257301434740       4.366263952193433E-003                        
   5.69397542244271       1.545391643363289E-003                        
   6.84397673183353       4.908361854610483E-004                        
   8.10392337664535       1.397086488908303E-004                        
   9.47531947589192       3.558361183644027E-005                        
   10.9598375637347       8.096466357554405E-006                        
   12.5593309474486       1.642708909566552E-006                        
   14.2758479323996       2.965938810009984E-007                        
   16.1116482030640       4.754732429514802E-008                        
   18.0692217104087       6.751169462255367E-009                        
   20.1513104920617       8.467187567022060E-010                        
   22.3609339469653       9.352020771970011E-011                        
   24.7014182063956       9.066606749167704E-012                        
   27.1764303961281       7.687363449753173E-013                        
   29.7900187807575       5.677534790338497E-014                        
   32.5466600353495       3.636338260986426E-015                        
   35.4513152220926       2.009810298315003E-016                        
   38.5094964891880       9.533662342373185E-018                        
   41.7273470970119       3.857758653674308E-019                        
   45.1117381723318       1.322590235296160E-020                        
   48.6703866831494       3.812520377375046E-022                        
   52.4120006467825       9.161202197099637E-024                        
   56.3464597342432       1.817193129409686E-025                        
   60.4850425305087       2.942483109948401E-027                        
   64.8407162572861       3.839941413589497E-029                        
   69.4285115890426       3.979093788538901E-031                        
   74.2660156887036       3.217746415876857E-033                        
   79.3740331847917       1.989360011844198E-035                        
   84.7774918932803       9.174803331715442E-038                        
   90.5067159166133       3.063597929422459E-040                        
   96.5992696619433       7.137877123560422E-043                        
   103.102726489260       1.107412963956401E-045                        
   110.079011673338       1.076644501393989E-048                        
   117.611597334412       6.044587863003915E-052                        
   125.818289122308       1.746746567183330E-055                        
   134.876188756134       2.186764947502585E-059                        
   145.077369431289       8.937404843361401E-064                        
   156.981623103640       6.961682775435475E-069                        
   172.032866745166       2.630674295402287E-075                        
*                                                                       
*  input file ao16p.scd                                                 
*                                                                       
                                                                        
1,16,                                                                   
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
                                                                        
0,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,4,                                                     
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   0.7708d0,                                                            
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
1,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,4,                                                     
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   0.6562d0,                                                            
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
2,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,3,                                                     
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   0.7708d0,                                                            
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
3,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,3,                                                     
   300.d0,                                                              
   300.d0,                                                              
   300.d0,                                                              
   0.6562d0,                                                            
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
4,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,2,                                                     
   300.d0,                                                              
   300.d0,                                                              
   0.7708d0,                                                            
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
5,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,2,                                                     
   300.d0,                                                              
   300.d0,                                                              
   0.6562d0,                                                            
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
6,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,1,                                                     
   300.d0,                                                              
   0.7708d0,                                                            
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
7,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,1,                                                     
   300.d0,                                                              
   0.6562d0,                                                            
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
8,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.7708d0,                                                            
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
9,                                                                      
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.6562d0,                                                            
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
10,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.4897d0,                                                            
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
11,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.3804d0,                                                            
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
12,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.2815d0,                                                            
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
13,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.21d0,                                                              
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
14,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.1541d0,                                                            
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
15,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.1129d0,                                                            
   0.0597d0,                                                            
   0.0312d0,                                                            
   0.0161d0,                                                            
   0.0083d0,                                                            
   0.0042d0,                                                            
   0.0021d0,                                                            
   0.0011d0,                                                            
   0.0005d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
16,                                                                     
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               
      -1.67709d0,  0.18963d0, 0.43546d0,                                
   0.318154709d0,0,                                                     
   0.0822d0,                                                            
   0.0432d0,                                                            
   0.0225d0,                                                            
   0.0116d0,                                                            
   0.0059d0,                                                            
   0.0030d0,                                                            
   0.0015d0,                                                            
   0.0008d0,                                                            
   0.0004d0,                                                            
   0.0002d0,                                                            
   0.0001d0,                                                            
   0.d0,                                                                
                                                                        
*                                                                       
*  output file ao.sr                                                    
*                                                                       
                                                                        
  Scattering States in Coulomb-like potential                           
                                                                        
                                                                        
  mass_1 =  4.02350000 mass_2 = 16.08800000                             
                                                                        
  Potential file name : ao16p.scd                                       
                                                                        
  lm =  3   omcst = 1.7500   b = 3.1396                                 
                                                                        
                                                                        
  Initial energy =      12.00                                           
                                                                        
                                                                        
  nmax=   24                                                            
                                                                        
  Coulomb phase shift eta_l   =    1.667727 rad                         
                      eta_l   =   95.553724 deg                         
                                                                        
  Colomb modified nuclear                                               
          phase shift delta_l =    1.508089 rad                         
                      delta_l =   86.407137 deg                         
                                                                        
  Coulomb modified nuclear                                              
  scattering amplitude   a'_l =    0.960157D-01 +i*  -0.727880D+00      
                                                                        
  wave function                                                         
                                                                        
    0.20 -0.2045D-03 -0.3257D-02                                        
    0.40 -0.2825D-02 -0.4500D-01                                        
    0.60 -0.1113D-01 -0.1773D+00                                        
    0.80 -0.2436D-01 -0.3879D+00                                        
    1.00 -0.3560D-01 -0.5670D+00                                        
    1.20 -0.3586D-01 -0.5711D+00                                        
    1.40 -0.2112D-01 -0.3364D+00                                        
    1.60  0.3525D-02  0.5615D-01                                        
    1.80  0.2618D-01  0.4170D+00                                        
    2.00  0.3527D-01  0.5617D+00                                        
    2.20  0.2631D-01  0.4190D+00                                        
    2.40  0.4077D-02  0.6493D-01                                        
    2.60 -0.2060D-01 -0.3282D+00                                        
    2.80 -0.3675D-01 -0.5853D+00                                        
    3.00 -0.3805D-01 -0.6060D+00                                        
    3.20 -0.2459D-01 -0.3916D+00                                        
    3.40 -0.1403D-02 -0.2234D-01                                        
    3.60  0.2445D-01  0.3894D+00                                        
    3.80  0.4671D-01  0.7439D+00                                        
    4.00  0.6146D-01  0.9789D+00                                        
    4.20  0.6743D-01  0.1074D+01                                        
    4.40  0.6530D-01  0.1040D+01                                        
    4.60  0.5683D-01  0.9051D+00                                        
    4.80  0.4402D-01  0.7011D+00                                        
    5.00  0.2869D-01  0.4569D+00                                        
                                                                        
                                                                        
