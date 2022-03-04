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
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
880D+00      
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
