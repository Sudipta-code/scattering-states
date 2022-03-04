ACHFCSPSE_SCATTERING.  CALCULATING SCATTERING STATES IN LOCAL AND       ACHF0000
1   NON-LOCAL COULOMB-LIKE POTENTIALS.  Z. PAPP.                        ACHF0000
REF. IN COMP. PHYS. COMMUN. 70 (1992) 435                               ACHF0000
*                                                                       ACHF0001
*  program cspse_scattering                                             ACHF0002
*     Zoltan Papp                                                       ACHF0003
*       present address (1991):                                         ACHF0004
*           Inst. Theor. Phys., Univ. Tubingen, Tubingen, Germany       ACHF0005
*       permanent address:                                              ACHF0006
*           Inst. Nucl. Res., Debrecen, Hungary                         ACHF0007
*=================================================================      ACHF0008
*                                                                       ACHF0009
*  potential separable expansion on Coulomb-Sturmian  basis             ACHF0010
*                                                                       ACHF0011
*  Calculating Scattering States in Local and Non-local                 ACHF0012
*  Coulomb-like Potentials                                              ACHF0013
*                                                                       ACHF0014
*=================================================================      ACHF0015
*                                                                       ACHF0016
        implicit real*8 (a-h,o-z)                                       ACHF0017
        character*12 twbp,glag00,glagm05                                ACHF0018
        parameter (ndim=25,npint=48)                                    ACHF0019
        common /iout/ iout,inp                                          ACHF0020
        common /intst/xst(npint),wst(npint),npst                        ACHF0021
        common /intho/xho(npint),who(npint),npho                        ACHF0022
        dimension vp(ndim,ndim)                                         ACHF0023
*                                                                       ACHF0024
        inp=5                                                           ACHF0025
        iout=6                                                          ACHF0026
        read(inp,*) glag00,glagm05                                      ACHF0027
        read(inp,*) xm1,xm2,twbp                                        ACHF0028
        read(inp,*) lm,omcst                                            ACHF0029
        read(inp,*) minn,maxn                                           ACHF0030
        if(maxn.ge.ndim) stop ' maxn > ndim  '                          ACHF0031
        if(minn.gt.maxn) stop ' minn > maxn '                           ACHF0032
*                                                                       ACHF0033
        open(2,file=glag00,status='old')                                ACHF0034
        read(2,*) npst,alfa                                             ACHF0035
        if(alfa.ne.0.0d0) stop ' alfa .ne. 0.0d0 '                      ACHF0036
        read(2,*) (xst(i),wst(i),i=1,npst)                              ACHF0037
        close(2)                                                        ACHF0038
*                                                                       ACHF0039
        open(2,file=glagm05,status='old')                               ACHF0040
        read(2,*) npho,alfa                                             ACHF0041
        if(alfa.ne.-0.5d0) stop ' alfa .ne. -0.5d0 '                    ACHF0042
        read(2,*) (xho(i),who(i),i=1,npho)                              ACHF0043
        close(2)                                                        ACHF0044
*                                                                       ACHF0045
        write(iout,'(/a,a/)') '  Scattering States',                    ACHF0046
     &             ' in Coulomb-like potential'                         ACHF0047
*                                                                       ACHF0048
        xm=xm1*xm2/(xm1+xm2)                                            ACHF0049
        write(iout,'(/1x,a,f12.8,a,f12.8//a,a)') ' mass_1 =',xm1,       ACHF0050
     &    ' mass_2 =',xm2,'  Potential file name : ',twbp               ACHF0051
        write(iout,'(/1x,a,i3,a,f7.4,a,f7.4)') ' lm =',lm,              ACHF0052
     &    '   omcst =',omcst,'   b =',sqrt(xm)*omcst                    ACHF0053
        call cspom(xm,lm,maxn,omcst,twbp,vp,z12,ipflag,1)               ACHF0054
        if(ipflag.lt.0) then                                            ACHF0055
          write(iout,'(a,i2)') ' The potential vanishes for l =',lm     ACHF0056
          stop                                                          ACHF0057
        endif                                                           ACHF0058
        call cstm(xm,lm,z12,omcst,vp,minn,maxn)                         ACHF0059
        stop                                                            ACHF0060
        end                                                             ACHF0061
                                                                        ACHF0062
        subroutine cstm(xm,lm,z12,omcst,vp,minn,maxn)                   ACHF0063
*                                                                       ACHF0064
*                                                                       ACHF0065
*  Formal parameters                                                    ACHF0066
*     xm     - reduced mass                                             ACHF0067
*     lm     - angular momentum                                         ACHF0068
*     omcst  - scalig parameter for the basis                           ACHF0069
*              b=sqrt(xm)*omcst                                         ACHF0070
*     z12    - charge                                                   ACHF0071
*     vp    -  matrix elements of the potential                         ACHF0072
*     minn  -  starting value of the radial quantum number              ACHF0073
*     maxn  -  maximum value of the radial quantum number               ACHF0074
*                                                                       ACHF0075
        implicit real*8 (a-h,o-z)                                       ACHF0076
        complex*16 vv,g,bwf,tl,fl,work,vdet                             ACHF0077
        parameter (ndim=25,hb2=41.801614d0,coul=1.43996518d0,           ACHF0078
     &          pi=3.141592653589793d0)                                 ACHF0079
        common /iout/ iout,inp                                          ACHF0080
        dimension szi(ndim),vp(ndim,ndim),vv(ndim,ndim),g(ndim,ndim),   ACHF0081
     &          bwf(ndim),cst(ndim),kpvt(ndim),work(ndim),vdet(2)       ACHF0082
        szfg(n,i)=(1-exp(-((i-n-1)*6.d0/n)**2))/(1-exp(-6.d0**2))       ACHF0083
*                                                                       ACHF0084
        bst=sqrt(xm)*omcst                                              ACHF0085
*                                                                       ACHF0086
*  <n'l| V_sl|nl>=<n'l| V |nl>-<n'l| V_c |nl>                           ACHF0087
*                                                                       ACHF0088
        do 1 i=1,maxn+1                                                 ACHF0089
           vp(i,i)=vp(i,i)-z12*coul                                     ACHF0090
   1    continue                                                        ACHF0091
*                                                                       ACHF0092
   5    read(inp,*,end=50) en,idr,dr                                    ACHF0093
          write(iout,'(//a,f10.2,/)') '  Initial energy = ',en          ACHF0094
*                                                                       ACHF0095
          do 20 nmm=minn,maxn                                           ACHF0096
            nmax=nmm                                                    ACHF0097
            ind=nmax+1                                                  ACHF0098
*                                                                       ACHF0099
*      calculation of the smoothing coefficients                        ACHF0100
*                                                                       ACHF0101
            do 12 i=1,ind                                               ACHF0102
  12          szi(i)=szfg(ind,i)                                        ACHF0103
            do 13 i=1,ind                                               ACHF0104
              do 13 k=1,i                                               ACHF0105
                 vv(k,i)=szi(i)*szi(k)*vp(k,i)                          ACHF0106
  13        continue                                                    ACHF0107
*                                                                       ACHF0108
*      V^{-1} - G(z)                                                    ACHF0109
*                                                                       ACHF0110
            call zsifa(vv,ndim,ind,kpvt,info)                           ACHF0111
            if(info.ne.0) stop ' singular potential matrix '            ACHF0112
            call zsidi(vv,ndim,ind,kpvt,vdet,work,1)                    ACHF0113
*                                                                       ACHF0114
            call csgrs(en,bst,xm,lm,z12,nmax,g,eta,cst)                 ACHF0115
            do 14 i=1,ind                                               ACHF0116
              do 14 j=1,i                                               ACHF0117
                 vv(j,i)=vv(j,i)-g(j,i)                                 ACHF0118
  14        continue                                                    ACHF0119
*                                                                       ACHF0120
*        cst( ) * ( V^{-1} - G(z) )^{-1} * cst( )                       ACHF0121
*                                                                       ACHF0122
            call zsifa(vv,ndim,ind,kpvt,info)                           ACHF0123
            if(info.ne.0) stop '  singular pse matrix '                 ACHF0124
            do 15 i=1,ind                                               ACHF0125
  15          bwf(i)=cst(i)                                             ACHF0126
            call zsisl(vv,ndim,ind,kpvt,bwf)                            ACHF0127
*                                                                       ACHF0128
            tl=(0.d0,0.d0)                                              ACHF0129
            do 16 i=1,ind                                               ACHF0130
  16          tl=tl+bwf(i)*cst(i)                                       ACHF0131
            zk=sqrt(2*xm/hb2*en)                                        ACHF0132
*                                                                       ACHF0133
* for complex potentials delta should be declaired as complex           ACHF0134
*                                                                       ACHF0135
            delta=log(1-2*(0.d0,1.d0)*zk*tl/en)/(2*(0.d0,1.d0))         ACHF0136
            write(iout,'(/a,i5)') '  nmax=',nmax                        ACHF0137
            write(iout,'(/a,f12.6,a/a,f12.6,a)')                        ACHF0138
     &          '  Coulomb phase shift eta_l   =',eta,' rad ',          ACHF0139
     &          '                      eta_l   =',eta*180/pi,' deg '    ACHF0140
            write(iout,'(/a/8x,a,f12.6,a/8x,a,f12.6,a)')                ACHF0141
     &          '  Colomb modified nuclear ',                           ACHF0142
     &          '  phase shift delta_l =',delta,' rad ',                ACHF0143
     &          '              delta_l =',delta*180/pi,' deg '          ACHF0144
            fl=-tl/en*exp(2*(0.d0,1.d0)*eta)                            ACHF0145
            write(iout,'(/a/a,d15.6,'' +i*'',d15.6)')                   ACHF0146
     &          '  Coulomb modified nuclear',                           ACHF0147
     &          '  scattering amplitude   a''_l = ',fl                  ACHF0148
  20      continue                                                      ACHF0149
*                                                                       ACHF0150
*       calculation of the wave function                                ACHF0151
*                                                                       ACHF0152
          if(idr.ne.0) then                                             ACHF0153
            write(iout,'(/a/)') '  wave function '                      ACHF0154
            call cswfs(en,bst,xm,lm,z12,nmax,bwf,idr,dr,1)              ACHF0155
          endif                                                         ACHF0156
*                                                                       ACHF0157
        goto 5                                                          ACHF0158
  50    return                                                          ACHF0159
        end                                                             ACHF0160
                                                                        ACHF0161
        subroutine cspom(xm,lm,nvmax,omcst,twbp,vp,z12,ipflag,job)      ACHF0162
*                                                                       ACHF0163
*   Calculates the Coulomb-Sturmian matrix elements of the total        ACHF0164
*   V = V_c + V_s  potential.                                           ACHF0165
*                                                                       ACHF0166
*  We suppose that the potential is a fish-bone optical potential       ACHF0167
*  which local part has the form                                        ACHF0168
*                                                                       ACHF0169
*  V_{loc}(r)=v0*(1+v1*r^2)*exp(-beta1*r^2)+v2*exp(-beta2*r^2)          ACHF0170
*      +zch*e^2*erf(beta3*r)/r,                                         ACHF0171
*                                                                       ACHF0172
*  and the nonlocal part is expressed on harmonic oscillator basis.     ACHF0173
*  We suppose also that the potential for high angular momentum         ACHF0174
*  becomes an angular momentum idependent local potential U_{loc}       ACHF0175
*                                                                       ACHF0176
*  The data file 'twbp' is organized as follows.                        ACHF0177
*  iswloc,lkmax  - control parameters: if(iswloc.eq.0) U_{loc}= 0,      ACHF0178
*         maximal angular momentum where V differ from U_{loc}          ACHF0179
*    zchr,v0r,beta1r,v1r,v2r,beta2r,beta3r - parameters for U_{loc}     ACHF0180
*  lk1     - angular momentum                                           ACHF0181
*    zch,v0,beta1,v1,v2,beta2,beta3 - parameters for the local part     ACHF0182
*              of the fish-bone potential                               ACHF0183
*    omfish,npfs - omega parameter of the harmonic oscillator basis     ACHF0184
*      and number of Pauli forbidden states in the fish-bone potential  ACHF0185
*    eta(i) - the first npfs eta(i) stand for the epsilon parameters,   ACHF0186
*             then eta(i) is the eta parameter. if(eta(i).eq.0) then    ACHF0187
*             there are no more non-zero eta parameters.                ACHF0188
*                                                                       ACHF0189
* Formal parameters:                                                    ACHF0190
*   xm     - reduced mass                                               ACHF0191
*   lm     - angular momentum                                           ACHF0192
*   nvmax  - maximal radial quantum number                              ACHF0193
*   omcst  - scalig parameter for the basis                             ACHF0194
*            b=sqrt(xm)*omcst                                           ACHF0195
*   twbp   - data file name for potential parameters                    ACHF0196
*   vp     - potential matrix                                           ACHF0197
*   z12    - charge                                                     ACHF0198
*   ipflag - control parameter                                          ACHF0199
*            ipflag=-2    V vanishes for lm                             ACHF0200
*            ipflag=-1    V vanishes for lm                             ACHF0201
*            ipflag= 0    V angular momentum independent                ACHF0202
*            ipflag= 1    V depends on the angular momentum lm          ACHF0203
*   job    - control parameter                                          ACHF0204
*            job=0    V-U_{loc} is calculated                           ACHF0205
*            job=1    V is calculated                                   ACHF0206
*                                                                       ACHF0207
        implicit real*8 (a-h,o-z)                                       ACHF0208
        real*8 mbar                                                     ACHF0209
        character*12 twbp                                               ACHF0210
        parameter (ndim=25,npint=48)                                    ACHF0211
        parameter (hb2=41.801614d0)                                     ACHF0212
*        parameter (hb2=1.0d0)                                          ACHF0213
        common /potpar/ zch,v0,v1,beta1,v2,beta2,beta3                  ACHF0214
        common /potpr/ zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r,iswloc     ACHF0215
        common /intho/ xho(npint),who(npint),npho                       ACHF0216
        common /intst/ xst(npint),wst(npint),npst                       ACHF0217
        dimension ho(ndim,npint),st(ndim,npint),eta(ndim),              ACHF0218
     &    mbar(ndim,ndim),tkii(ndim),tkmi(ndim),tr(ndim,ndim),          ACHF0219
     &    vfb(ndim,ndim),vp(ndim,ndim),pf(npint)                        ACHF0220
*                                                                       ACHF0221
* reading in the potential parameters                                   ACHF0222
*                                                                       ACHF0223
        open(2,file=twbp,status='old')                                  ACHF0224
        read(2, *) iswloc,lkmax                                         ACHF0225
        read(2, *) zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r                ACHF0226
*                                                                       ACHF0227
   1    read(2, *,end=4) lk1                                            ACHF0228
        read(2, *) zch,v0,v1,beta1,v2,beta2,beta3                       ACHF0229
*                                                                       ACHF0230
        read(2,*) omfish,npfs                                           ACHF0231
        i=0                                                             ACHF0232
   2    i=i+1                                                           ACHF0233
        read(2,*) eta(i)                                                ACHF0234
        if(i.le.npfs.or.eta(i).ne.0.d0) goto 2                          ACHF0235
        neta=i-1                                                        ACHF0236
        if(lk1.ne.lm) then                                              ACHF0237
          goto 1                                                        ACHF0238
        else                                                            ACHF0239
          close(2)                                                      ACHF0240
          ipflag=2                                                      ACHF0241
          z12=zch                                                       ACHF0242
          goto 5                                                        ACHF0243
        endif                                                           ACHF0244
   4    close(2)                                                        ACHF0245
        if(job.eq.1) then                                               ACHF0246
          if(iswloc.ne.0) then                                          ACHF0247
            ipflag=1                                                    ACHF0248
            z12=zchr                                                    ACHF0249
            neta=0                                                      ACHF0250
          else                                                          ACHF0251
            ipflag=0                                                    ACHF0252
            return                                                      ACHF0253
          endif                                                         ACHF0254
        elseif(job.eq.0) then                                           ACHF0255
          if(lm.lt.lkmax) then                                          ACHF0256
            ipflag=-1                                                   ACHF0257
          else                                                          ACHF0258
            ipflag=-2                                                   ACHF0259
          endif                                                         ACHF0260
          return                                                        ACHF0261
        else                                                            ACHF0262
          stop  ' error in calling cspom '                              ACHF0263
        endif                                                           ACHF0264
   5    continue                                                        ACHF0265
*                                                                       ACHF0266
*  calculation of the matrix elements of the local part of the          ACHF0267
*  potential on Coulomb-Sturmian basis                                  ACHF0268
*                                                                       ACHF0269
        nn=nvmax+1                                                      ACHF0270
        bst=sqrt(xm)*omcst                                              ACHF0271
        do 6 i=1,npst                                                   ACHF0272
          call csturm(st(1,i),xst(i),lm,nvmax)                          ACHF0273
          xib=xst(i)/(2*bst)                                            ACHF0274
          if(job.eq.1) then                                             ACHF0275
            if(ipflag.eq.1) then                                        ACHF0276
              pf(i)=potfr(xib)/(2*bst)                                  ACHF0277
            else                                                        ACHF0278
              pf(i)=potfun(xib)/(2*bst)                                 ACHF0279
            endif                                                       ACHF0280
          else                                                          ACHF0281
            if(iswloc.eq.0) then                                        ACHF0282
              pf(i)=potfun(xib)/(2*bst)                                 ACHF0283
            else                                                        ACHF0284
              pf(i)=(potfun(xib)-potfr(xib))/(2*bst)                    ACHF0285
            endif                                                       ACHF0286
          endif                                                         ACHF0287
   6    continue                                                        ACHF0288
        do 8 i=1,nn                                                     ACHF0289
          do 8 j=1,i                                                    ACHF0290
            vv=0.0d0                                                    ACHF0291
            do 7 n=1,npst                                               ACHF0292
   7           vv=vv+st(i,n)*st(j,n)*pf(n)*wst(n)                       ACHF0293
            vp(j,i)=vv                                                  ACHF0294
            vp(i,j)=vv                                                  ACHF0295
   8    continue                                                        ACHF0296
*                                                                       ACHF0297
        if(neta.ne.0) then                                              ACHF0298
*                                                                       ACHF0299
*  calculation of the matrix elements of the non-local part of the      ACHF0300
*  potential on Coulomb-Sturmian basis                                  ACHF0301
*                                                                       ACHF0302
*   calculate the norm kernel eigenfunctions                            ACHF0303
*   (harmonic oscillator functions in tis case)                         ACHF0304
*                                                                       ACHF0305
          bho=1.d0/sqrt(xm*omfish)                                      ACHF0306
          do 9 i=1,npho                                                 ACHF0307
            xsq=sqrt(xho(i))                                            ACHF0308
            call hof(ho(1,i),xsq,lm,nvmax)                              ACHF0309
            pf(i)=potfun(xsq*bho)                                       ACHF0310
   9      continue                                                      ACHF0311
          do 10 i=neta+1,nn                                             ACHF0312
  10        eta(i)=0.d0                                                 ACHF0313
*                                                                       ACHF0314
*   calculate the M-bar matrix and the matrix elements of the           ACHF0315
*   kinetic energy                                                      ACHF0316
*                                                                       ACHF0317
          hbm=hb2/2*omfish                                              ACHF0318
          do 12 j=1,nn                                                  ACHF0319
          do 12 i=1,j                                                   ACHF0320
          if(i.le.npfs) then                                            ACHF0321
             mbar(i,j)=1.d0                                             ACHF0322
          else                                                          ACHF0323
            if(i.eq.j.or.i.gt.neta) then                                ACHF0324
             mbar(i,j)=0.d0                                             ACHF0325
            else                                                        ACHF0326
             mbar(i,j)=1.d0-sqrt((1.d0-eta(i))/(1.d0-eta(j)))           ACHF0327
            endif                                                       ACHF0328
          endif                                                         ACHF0329
          if(i.le.neta) then                                            ACHF0330
            if(i.eq.j) then                                             ACHF0331
              tkii(i)=hbm*(2*i+lm-0.5d0)                                ACHF0332
              if(i.le.npfs) tkii(i)=tkii(i)-eta(i)                      ACHF0333
            elseif(i.eq.j-1) then                                       ACHF0334
              tkmi(i)=hbm*sqrt(i*(i+lm+0.5d0))                          ACHF0335
            endif                                                       ACHF0336
          endif                                                         ACHF0337
  12      continue                                                      ACHF0338
*                                                                       ACHF0339
*   calcullate the harmonic oscillator matrix elements of the           ACHF0340
*   non-local part of the potential                                     ACHF0341
*                                                                       ACHF0342
          do 17 j=1,nn                                                  ACHF0343
          do 17 i=1,j                                                   ACHF0344
          if(mbar(i,j).eq.0.d0) then                                    ACHF0345
            vfb(i,j)=0.0d0                                              ACHF0346
          else                                                          ACHF0347
            vv=0.0d0                                                    ACHF0348
            do 15 n=1,npho                                              ACHF0349
  15        vv=vv+ho(i,n)*ho(j,n)*pf(n)*who(n)                          ACHF0350
            if(i.eq.j) then                                             ACHF0351
              vv=vv+tkii(i)                                             ACHF0352
            elseif(i.eq.j-1) then                                       ACHF0353
              vv=vv+tkmi(i)                                             ACHF0354
            endif                                                       ACHF0355
            vfb(i,j)=vv*mbar(i,j)                                       ACHF0356
          endif                                                         ACHF0357
          vfb(j,i)=vfb(i,j)                                             ACHF0358
  17      continue                                                      ACHF0359
*                                                                       ACHF0360
*  transform the matrix elemnts from harmonic oscillator basis          ACHF0361
*  to Coulomb-Sturmian basis                                            ACHF0362
*                                                                       ACHF0363
          call csho(lm,nvmax,bst,bho,tr)                                ACHF0364
*                                                                       ACHF0365
          do 20 j=1,nn                                                  ACHF0366
          do 20 i=1,j                                                   ACHF0367
            vv=0.0d0                                                    ACHF0368
            do 19 k1=1,nn                                               ACHF0369
            do 19 k2=1,nn                                               ACHF0370
               vv=vv+tr(i,k1)*vfb(k1,k2)*tr(j,k2)                       ACHF0371
  19        continue                                                    ACHF0372
            vp(i,j)=vp(i,j)-vv                                          ACHF0373
            vp(j,i)=vp(i,j)                                             ACHF0374
  20      continue                                                      ACHF0375
        endif                                                           ACHF0376
*                                                                       ACHF0377
        return                                                          ACHF0378
        end                                                             ACHF0379
                                                                        ACHF0380
        subroutine csho(l,nvmax,bst,bho,tr)                             ACHF0381
*                                                                       ACHF0382
*   Calculates the overlap matrix elements betwen Coulomb-Sturmian      ACHF0383
*   and norm kernel eigenfunctions (harmonic  oscillator functions)     ACHF0384
*                                                                       ACHF0385
*   l     - angular momentum                                            ACHF0386
*   nvmax - maximal radial quantum number                               ACHF0387
*   bst   - b parameter of the Coulomb-Sturmian functions               ACHF0388
*   bho   - b parameter of the harmonic oscillator functions            ACHF0389
*   tr    - overlap metrix elements                                     ACHF0390
*                                                                       ACHF0391
        implicit real*8 (a-h,o-z)                                       ACHF0392
        parameter (ndim=25,npint=48)                                    ACHF0393
        common /intst/xst(npint),wst(npint),npst                        ACHF0394
        dimension  h(ndim,npint),s(ndim,npint),tr(ndim,ndim)            ACHF0395
        do 1 i=1,npst                                                   ACHF0396
        call csturm(s(1,i),2*xst(i),l,nvmax)                            ACHF0397
        xib=xst(i)/bst                                                  ACHF0398
        call howf(h(1,i),xib,bho,l,nvmax)                               ACHF0399
   1    continue                                                        ACHF0400
        do 8 ni=1,nvmax+1                                               ACHF0401
        do 8 nj=1,nvmax+1                                               ACHF0402
        vv=0.d0                                                         ACHF0403
        do 3 i=1,npst                                                   ACHF0404
   3    vv=vv+wst(i)*s(ni,i)*h(nj,i)                                    ACHF0405
        tr(ni,nj)=vv/bst                                                ACHF0406
   8    continue                                                        ACHF0407
        return                                                          ACHF0408
        end                                                             ACHF0409
                                                                        ACHF0410
                                                                        ACHF0411
        subroutine csturm(st,x,l,nmax)                                  ACHF0412
*                                                                       ACHF0413
*   Calculates Coulomb-Sturmian fuctions                                ACHF0414
*   (without the exponential factor)                                    ACHF0415
*                                                                       ACHF0416
*   st(n+1)=[n!/(n+2l+1)!]^{1/2} x^{l+1} L_n^{2l+1} (x)                 ACHF0417
*           from n=0 to nmax                                            ACHF0418
*   x  - radial variable (2b*r)                                         ACHF0419
*   l  - angular momentum                                               ACHF0420
*                                                                       ACHF0421
        implicit real*8 (a-h,o-z)                                       ACHF0422
        dimension st(1)                                                 ACHF0423
        l2=2*l+1                                                        ACHF0424
        sfl=1.d0                                                        ACHF0425
        do 1 i=2,l2                                                     ACHF0426
   1    sfl=sfl*i                                                       ACHF0427
        xl0=0.d0                                                        ACHF0428
        st(1)=x**(l+1)/sqrt(sfl)                                        ACHF0429
        do 5 i=1,nmax                                                   ACHF0430
        df=sqrt(dble(i*(i+l2)))                                         ACHF0431
        st(i+1)=((2*i+l2-1-x)*st(i)-xl0)/df                             ACHF0432
        xl0=df*st(i)                                                    ACHF0433
   5    continue                                                        ACHF0434
        return                                                          ACHF0435
        end                                                             ACHF0436
                                                                        ACHF0437
        subroutine hof(ho,x,l,nmax)                                     ACHF0438
*                                                                       ACHF0439
*  ho(n+1)=howf(n,l,r,r0)*sqrt(r0/2)*exp((r/r0)**2/2)                   ACHF0440
*  where howf is the harmonic oscillator function with                  ACHF0441
*  parameters                                                           ACHF0442
*   n:  radial quantum number                                           ACHF0443
*   l:  angular momentum variable                                       ACHF0444
*   r:  radial variable                                                 ACHF0445
*   r0: oscillator lenght parameter                                     ACHF0446
*  --Other input parameters--                                           ACHF0447
*   x:  r/r0                                                            ACHF0448
*   nn: maximal radial quantum number                                   ACHF0449
*  --Output parameters--                                                ACHF0450
*   ho: howf(n,l,r,r0)*sqrt(r0/2)*exp((r/r0)**2/2)                      ACHF0451
*                                                                       ACHF0452
        implicit real*8 (a-h,o-z)                                       ACHF0453
        parameter (sqpi=1.772453850905516d0)                            ACHF0454
        dimension ho(1)                                                 ACHF0455
        x2=x*x                                                          ACHF0456
        cl=2*x2/sqpi                                                    ACHF0457
        do 1 i=1,l                                                      ACHF0458
    1   cl=2*cl/(2*i+1)*x2                                              ACHF0459
        ho(1)=sqrt(cl)                                                  ACHF0460
        rn1=0.d0                                                        ACHF0461
        do 2 i=1,nmax                                                   ACHF0462
        df=sqrt(i*(i+l+0.5d0))                                          ACHF0463
        ho(i+1)=((2*i+l-0.5d0-x2)*ho(i)-rn1)/df                         ACHF0464
        rn1=df*ho(i)                                                    ACHF0465
    2   continue                                                        ACHF0466
        return                                                          ACHF0467
        end                                                             ACHF0468
                                                                        ACHF0469
        subroutine howf(ho,r,r0,l,nmax)                                 ACHF0470
*                                                                       ACHF0471
*  ho(n+1)=howf(n,l,r,r0)                                               ACHF0472
*  where howf is the harmonic oscillator function with                  ACHF0473
*  parameters                                                           ACHF0474
*    n:  radial quantum number                                          ACHF0475
*    l:  angular momentum variable                                      ACHF0476
*    r:  radial variable                                                ACHF0477
*    r0: oscillator lenght parameter                                    ACHF0478
*  --Other input parameters--                                           ACHF0479
*  nn: maximal radial quantum number                                    ACHF0480
*  --Output parameters--                                                ACHF0481
*  ho: howf(n,l,r,r0)                                                   ACHF0482
*                                                                       ACHF0483
        implicit real*8 (a-h,o-z)                                       ACHF0484
        dimension ho(1)                                                 ACHF0485
        x=r/r0                                                          ACHF0486
        aa=sqrt(2/r0)*exp(-x**2/2)                                      ACHF0487
        call hof(ho,x,l,nmax)                                           ACHF0488
        do 1 i=0,nmax                                                   ACHF0489
   1    ho(i+1)=aa*ho(i+1)                                              ACHF0490
        return                                                          ACHF0491
        end                                                             ACHF0492
                                                                        ACHF0493
        function potfun(r)                                              ACHF0494
*                                                                       ACHF0495
*  The potential parameters are stoored in the commonn block /potpar/.  ACHF0496
*  In this case they describe the following potential:                  ACHF0497
*  rpot=v0*(1+v1*r^2)*exp(-beta1*r^2)+v2*exp(-beta2*r^2)                ACHF0498
*      +zch*coul*erf(beta3*r)/r                                         ACHF0499
*                                                                       ACHF0500
        implicit real*8 (a-h,o-z)                                       ACHF0501
        common /potpar/ zch,v0,v1,beta1,v2,beta2,beta3                  ACHF0502
        common /potpr/ zchr,v0r,v1r,beta1r,v2r,beta2r,beta3r,iswloc     ACHF0503
        parameter (coul=1.43996518d0)                                   ACHF0504
*        parameter (coul=1.0d0)                                         ACHF0505
*                                                                       ACHF0506
        rr=r*r                                                          ACHF0507
        rpot=v0*(1.d0+v1*rr)*exp(-beta1*rr)                             ACHF0508
     &          +v2*exp(-beta2*rr)                                      ACHF0509
        if(zch.eq.0.d0) then                                            ACHF0510
          rpotc=0.d0                                                    ACHF0511
        else                                                            ACHF0512
          rpotc=zch*coul*erf(beta3*r)/r                                 ACHF0513
        endif                                                           ACHF0514
        potfun=rpot+rpotc                                               ACHF0515
        return                                                          ACHF0516
*                                                                       ACHF0517
        entry potfr(r)                                                  ACHF0518
*                                                                       ACHF0519
*       local auxiliary potential U_{loc}                               ACHF0520
*                                                                       ACHF0521
*  The potential parameters are stoored in the commonn block /potpr/.   ACHF0522
*  In this case they describe the following potential:                  ACHF0523
*  rpot=v0r*(1+v1r*r^2)*exp(-beta1r*r^2)+v2r*exp(-beta2r*r^2)           ACHF0524
*      +zchr*coul*erf(beta3r*r)/r                                       ACHF0525
*                                                                       ACHF0526
        rr=r*r                                                          ACHF0527
        rpot=v0r*(1.d0+v1r*rr)*exp(-beta1r*rr)                          ACHF0528
     &          +v2r*exp(-beta2r*rr)                                    ACHF0529
        if(zchr.eq.0.d0) then                                           ACHF0530
          rpotc=0.d0                                                    ACHF0531
        else                                                            ACHF0532
          rpotc=zchr*coul*erf(beta3r*r)/r                               ACHF0533
        endif                                                           ACHF0534
        potfr=rpot+rpotc                                                ACHF0535
        return                                                          ACHF0536
        end                                                             ACHF0537
                                                                        ACHF0538
        function erf(x)                                                 ACHF0539
*                                                                       ACHF0540
*       calculates the error function  (accuracy > 1.2d-7 )             ACHF0541
*                                                                       ACHF0542
        implicit real*8 (a-h,o-z)                                       ACHF0543
        z=abs(x)                                                        ACHF0544
        t=1.d0/(1.d0+0.5d0*z)                                           ACHF0545
        erfc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+                   ACHF0546
     &    t*(0.37409196d0+t*(0.09678418d0+t*(-0.18628806d0+             ACHF0547
     &    t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0+             ACHF0548
     &    t*(-0.82215223d0+t*0.17087277)))))))))                        ACHF0549
        if(x.lt.0.d0) erfc=2.d0-erfc                                    ACHF0550
        erf=1.d0-erfc                                                   ACHF0551
        return                                                          ACHF0552
        end                                                             ACHF0553
                                                                        ACHF0554
        subroutine csgrs(e,bst,xm,l,z12,nmax,g,eta,cst)                 ACHF0555
*                                                                       ACHF0556
*       Calculates the matrix elements of the                           ACHF0557
*       Coulomb-Green operator on Coulomb-Sturmian basis                ACHF0558
*                                                                       ACHF0559
*     e      - energy parameter of the Green operator                   ACHF0560
*     bst    - b parameter for the basis                                ACHF0561
*     xm     - reduced mass                                             ACHF0562
*     lm     - angular momentum                                         ACHF0563
*     z12    - charge                                                   ACHF0564
*     nmax   - maximum value of the radial quantum number               ACHF0565
*     g      - matrix elements of the Green operator                    ACHF0566
*                                                                       ACHF0567
        implicit real*8 (a-h,o-z)                                       ACHF0568
        complex*16 g,zki,gammai,s0,f21,clogam                           ACHF0569
        parameter (ndim=25,hb2=41.801614d0,coul=1.43996518d0,           ACHF0570
     &          pi=3.141592653589793d0)                                 ACHF0571
        dimension g(ndim,ndim),cst(ndim)                                ACHF0572
*                                                                       ACHF0573
        zk=sqrt(2*xm/hb2*e)                                             ACHF0574
        zki=(0.0d0,1.0d0)*zk                                            ACHF0575
        ctr=z12*coul*xm/hb2                                             ACHF0576
        gamma=ctr/zk                                                    ACHF0577
        gammai=(0.0d0,1.0d0)*gamma                                      ACHF0578
*                                                                       ACHF0579
*  < 0l | \Delta g_l^c (E) \Delta | 0l>                                 ACHF0580
*                                                                       ACHF0581
        g(1,1)=-4*xm*bst/(hb2*(l+gammai+1)*(bst-zki)**2)                ACHF0582
     &          *f21(-l+gammai,l+gammai+2,((bst+zki)/(bst-zki))**2)     ACHF0583
*                                                                       ACHF0584
        zkb2=zk**2+bst**2                                               ACHF0585
        sx1=(zk**2-bst**2)/zkb2                                         ACHF0586
        sx2=4*ctr*bst/zkb2                                              ACHF0587
        s0=4*xm/hb2*bst/zkb2                                            ACHF0588
*                                                                       ACHF0589
*  < 0l | \Delta g_l^c (E) \Delta | il>                                 ACHF0590
*                                                                       ACHF0591
        do 5 i=0,nmax-1                                                 ACHF0592
          df=sqrt(dble((i+1)*(i+2*l+2)))                                ACHF0593
          g(1,i+2)=(g(1,i+1)*(2*(i+l+1)*sx1-sx2)-s0)/df                 ACHF0594
          g(i+2,1)=g(1,i+2)                                             ACHF0595
          s0=g(1,i+1)*df                                                ACHF0596
   5    continue                                                        ACHF0597
*                                                                       ACHF0598
*  < jl | \Delta g_l^c (E) \Delta | il>                                 ACHF0599
*                                                                       ACHF0600
        do 6 i=0,nmax-1                                                 ACHF0601
          s0=(0.d0,0.d0)                                                ACHF0602
          do 6 j=0,i                                                    ACHF0603
            df=sqrt(dble((j+1)*(j+2*l+2)))                              ACHF0604
            g(j+2,i+2)=(g(j+1,i+2)*(2*(j+l+1)*sx1-sx2)-s0)/df           ACHF0605
            g(i+2,j+2)=g(j+2,i+2)                                       ACHF0606
            s0=g(j+1,i+2)*df                                            ACHF0607
   6    continue                                                        ACHF0608
*                                                                       ACHF0609
        zb=zk/bst                                                       ACHF0610
        if(z12.eq.0.0d0) then                                           ACHF0611
          eta=0.d0                                                      ACHF0612
          cst(1)=(2*zb/(1+zb*zb))**(l+1)                                ACHF0613
        else                                                            ACHF0614
          s0=clogam(l+gammai+1)                                         ACHF0615
          eta=dimag(s0)                                                 ACHF0616
          cst(1)=sqrt(2*pi*gamma/(exp(2*pi*gamma)-1))*                  ACHF0617
     &           exp(2*gamma*atan(zb))*(2*zb/(1+zb*zb))**(l+1)          ACHF0618
        endif                                                           ACHF0619
        do 7 i=1,l                                                      ACHF0620
   7      cst(1)=cst(1)*sqrt((gamma*gamma+i*i)/(i*(i+0.5d0)))           ACHF0621
        s2=0.d0                                                         ACHF0622
        do 8 i=0,nmax-1                                                 ACHF0623
          df=sqrt(dble((i+1)*(i+2*l+2)))                                ACHF0624
          cst(i+2)=(cst(i+1)*(2*(i+l+1)*sx1-sx2)-s2)/df                 ACHF0625
          s2=cst(i+1)*df                                                ACHF0626
   8    continue                                                        ACHF0627
*                                                                       ACHF0628
        return                                                          ACHF0629
        end                                                             ACHF0630
                                                                        ACHF0631
        subroutine cswfs(e,bst,xm,l,z12,nmax,bwf,idr,dr,isw)            ACHF0632
*                                                                       ACHF0633
*       Calculates the wave function                                    ACHF0634
*                                                                       ACHF0635
*     e      - energy parameter of the Coulomb-Green operator           ACHF0636
*     bst    - b parameter for the basis                                ACHF0637
*     xm     - reduced mass                                             ACHF0638
*     lm     - angular momentum                                         ACHF0639
*     z12    - charge                                                   ACHF0640
*     nmax   - maximum value of the radial quantum number               ACHF0641
*     bwf    - wave function coefficients b                             ACHF0642
*     idr    - number of points in calculating wave function            ACHF0643
*     dr     - step size in calculating wave function                   ACHF0644
*     isw    - if(isw.ne.0) the Coulomb function is added on            ACHF0645
*                                                                       ACHF0646
        implicit real*8 (a-h,o-z)                                       ACHF0647
        complex*16 bwf,zki,gammai,xum,cheby,gcs,gexpi,a0,s              ACHF0648
        parameter (ndim=25,ldim=100,eps=1.d-15,step=100)                ACHF0649
        parameter (hb2=41.801614d0,coul=1.43996518d0)                   ACHF0650
        common /iout/ iout                                              ACHF0651
        dimension bwf(ndim),cheby(ldim),gcs(ndim),                      ACHF0652
     &    fc(2*ndim),fcp(1),gc(1),gcp(1)                                ACHF0653
        equivalence (gcs(1),fc(1))                                      ACHF0654
*                                                                       ACHF0655
        sfl=1                                                           ACHF0656
        do 1 i=1,2*l+1                                                  ACHF0657
   1      sfl=sfl*i                                                     ACHF0658
        sfl=(2*bst)**(l+1)/sqrt(sfl)                                    ACHF0659
        zk=sqrt(2*xm/hb2*e)                                             ACHF0660
        zki=(0.d0,1.d0)*zk                                              ACHF0661
        ctr=z12*coul*xm/hb2                                             ACHF0662
        gamma=ctr/zk                                                    ACHF0663
        gammai=(0.0d0,1.0d0)*gamma                                      ACHF0664
*                                                                       ACHF0665
        zkb2=zk**2+bst**2                                               ACHF0666
        sx1=(zk**2-bst**2)/zkb2                                         ACHF0667
        sx2=4*ctr*bst/zkb2                                              ACHF0668
        s0=4*xm/hb2*bst/zkb2                                            ACHF0669
*                                                                       ACHF0670
        xum=-2*xm/hb2*sfl/(bst-zki)*                                    ACHF0671
     &          exp((l-gammai)*log(2*zki/(zki-bst)))                    ACHF0672
        call ch1f0(-l+gammai,-(bst+zki)/(2*zki),ldim,cheby)             ACHF0673
        nchu=10                                                         ACHF0674
        do 8 np=1,idr                                                   ACHF0675
          r=np*dr                                                       ACHF0676
          rho=zk*r                                                      ACHF0677
*                                                                       ACHF0678
*  < r |  g_l^c (E) \Delta | 0l>                                        ACHF0679
*                                                                       ACHF0680
          gcs(1)=xum*r**(l+1)*exp(-bst*r)*                              ACHF0681
     &          gexpi(-(bst+zki)*r,l+gammai,cheby,nchu)                 ACHF0682
          st0=0.d0                                                      ACHF0683
          df=0.d0                                                       ACHF0684
          st1=sfl*exp(-bst*r)*r**(l+1)                                  ACHF0685
          a0=(0.d0,0.d0)                                                ACHF0686
*                                                                       ACHF0687
*  < r |  g_l^c (E) \Delta | il>                                        ACHF0688
*                                                                       ACHF0689
          do 4 i=0,nmax-1                                               ACHF0690
            st2=2*(i+l+1-bst*r)*st1-df*st0                              ACHF0691
            df=sqrt(dble((i+1)*(i+2*l+2)))                              ACHF0692
            gcs(i+2)=(gcs(i+1)*(2*(i+l+1)*sx1-sx2)-a0-s0*st1)/df        ACHF0693
            a0=df*gcs(i+1)                                              ACHF0694
            st0=st1                                                     ACHF0695
            st1=st2/df                                                  ACHF0696
   4      continue                                                      ACHF0697
*                                                                       ACHF0698
*  < r| u_l>=\sum_{i=0}^{nmax} bwf(i+1) * <r|g_l^c (E) \Delta|il>       ACHF0699
*                                                                       ACHF0700
          s=(0.d0,0.d0)                                                 ACHF0701
          do 5 i=1,nmax+1                                               ACHF0702
   5        s=s+bwf(i)*gcs(i)                                           ACHF0703
*                                                                       ACHF0704
          if(isw.ne.0) then                                             ACHF0705
             call rcwff(rho,gamma,l,l,fc,gc,fcp,gcp,eps,step,1)         ACHF0706
             s=fc(l+1)+s                                                ACHF0707
          endif                                                         ACHF0708
*                                                                       ACHF0709
          write(iout,'(2x,f6.2,2d12.4)')  r,s                           ACHF0710
*                                                                       ACHF0711
   8    continue                                                        ACHF0712
*                                                                       ACHF0713
        return                                                          ACHF0714
        end                                                             ACHF0715
                                                                        ACHF0716
        function f21(a,c,z)                                             ACHF0717
*                                                                       ACHF0718
*  calculates the hypergeometric function 2_F_1 ( a, 1; c; z )          ACHF0719
*                                                                       ACHF0720
        complex*16 f21,a,c,z,p0,q0,q1,r,p2,q2,qq                        ACHF0721
        real*8 eps                                                      ACHF0722
        parameter (eps=1.d-14,nn=200)                                   ACHF0723
        common /iout/ iout                                              ACHF0724
*                                                                       ACHF0725
        p0=(1.d0,0.d0)                                                  ACHF0726
        q0=(0.d0,0.d0)                                                  ACHF0727
        q1=(1.d0,0.d0)                                                  ACHF0728
        qq=(0.d0,0.d0)                                                  ACHF0729
        n=0                                                             ACHF0730
   1    if(abs(1-qq/q1).gt.eps) then                                    ACHF0731
           r=-(a+n)*(c+n-1)/((c+2*n-1)*(c+2*n))*z                       ACHF0732
           p2=1+r*p0                                                    ACHF0733
           q2=q1+r*q0                                                   ACHF0734
           p0=1/p2                                                      ACHF0735
           q0=q1*p0                                                     ACHF0736
           qq=q2*p0                                                     ACHF0737
           n=n+1                                                        ACHF0738
           r=-n*(c+n-a-1)/((c+2*n-2)*(c+2*n-1))*z                       ACHF0739
           p2=1+r*p0                                                    ACHF0740
           q2=qq+r*q0                                                   ACHF0741
           p0=1/p2                                                      ACHF0742
           q0=qq*p0                                                     ACHF0743
           q1=q2*p0                                                     ACHF0744
           if(n.gt.nn) then                                             ACHF0745
              f21=q1                                                    ACHF0746
              write(iout,'('' estimated accuracy in f21 '',d16.9)')     ACHF0747
     &          abs((qq-q1)/q1)                                         ACHF0748
              return                                                    ACHF0749
           endif                                                        ACHF0750
           goto 1                                                       ACHF0751
        else                                                            ACHF0752
           f21=q1                                                       ACHF0753
           return                                                       ACHF0754
        endif                                                           ACHF0755
        end                                                             ACHF0756
                                                                        ACHF0757
        function clogam(z)                                              ACHF0758
c                                                                       ACHF0759
c    computes the logarithm of the gamma function                       ACHF0760
c    for any complex argument K.S Kolbig, CPC 4 (1972) 221              ACHF0761
c                                                                       ACHF0762
        implicit real*8 (a-h,o-z)                                       ACHF0763
        complex*16 z,v,h,r,clogam                                       ACHF0764
        common /iout/ iout                                              ACHF0765
        dimension b(10)                                                 ACHF0766
c                                                                       ACHF0767
        data lerr /3/                                                   ACHF0768
        data pi /3.14159 26535 89793d0/                                 ACHF0769
        data b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),b(9),b(10)         ACHF0770
     &       /+8.33333 33333 3333d-2, -2.77777 77777 7778d-3,           ACHF0771
     &        +7.93650 79365 0794d-4, -5.95238 09523 8095d-4,           ACHF0772
     &        +8.41750 84175 0842d-4, -1.91752 69175 2692d-3,           ACHF0773
     &        +6.41025 64102 5641d-3, -2.95506 53594 7712d-2,           ACHF0774
     &        +1.79644 37236 8831d-1, -1.39243 22169 0590d+0/           ACHF0775
c                                                                       ACHF0776
        x=dble(z)                                                       ACHF0777
        t=dimag(z)                                                      ACHF0778
        if(-abs(x) .eq. aint(x) .and. t .eq. 0.d0) go to 5              ACHF0779
        f=abs(t)                                                        ACHF0780
        v=dcmplx(x,f)                                                   ACHF0781
        if(x .lt. 0.d0) v=1.d0-v                                        ACHF0782
        h=0.d0                                                          ACHF0783
        c=dble(v)                                                       ACHF0784
        if(c .ge. 7.d0) go to 3                                         ACHF0785
        n=6-int(c)                                                      ACHF0786
        h=v                                                             ACHF0787
        d=dimag(v)                                                      ACHF0788
        a=atan2(d,c)                                                    ACHF0789
        if(n .eq. 0) go to 2                                            ACHF0790
        do 1 i = 1,n                                                    ACHF0791
        c=c+1.d0                                                        ACHF0792
        v=dcmplx(c,d)                                                   ACHF0793
        h=h*v                                                           ACHF0794
   1    a=a+atan2(d,c)                                                  ACHF0795
   2    h=dcmplx(0.5d0*log(dble(h)**2+dimag(h)**2),a)                   ACHF0796
        v=v+1.d0                                                        ACHF0797
   3    r=1.d0/v**2                                                     ACHF0798
        clogam=.918938533204673d0+(v-.5d0)*log(v)-v+(b(1)+r*(b(2)+      ACHF0799
     &  r*(b(3)+r*(b(4)+r*(b(5)+r*(b(6)+r*(b(7)+r*(b(8)+r*(b(9)+        ACHF0800
     &  r*b(10))))))))))/v-h                                            ACHF0801
        if(x .ge. 0.d0) go to 4                                         ACHF0802
c                                                                       ACHF0803
        a=aint(x)-1.d0                                                  ACHF0804
        c=pi*(x-a)                                                      ACHF0805
        d=pi*f                                                          ACHF0806
        e=exp(-2.d0*d)                                                  ACHF0807
        f=sin(c)                                                        ACHF0808
        e=d+0.5d0*log(e*f**2+0.25d0*(1.d0-e)**2)                        ACHF0809
        f=atan2(cos(c)*tanh(d),f)-a*pi                                  ACHF0810
        clogam=1.1447298858494d0-dcmplx(e,f)-clogam                     ACHF0811
c                                                                       ACHF0812
   4    if(sign(1.d0,t) .eq. -1.d0) clogam=conjg(clogam)                ACHF0813
        return                                                          ACHF0814
c                                                                       ACHF0815
   5    write(iout,100) x                                               ACHF0816
        clogam=(0.d0,0.d0)                                              ACHF0817
        return                                                          ACHF0818
 100    format(' clogam ... argument is non positive integer = ',f20.2) ACHF0819
c                                                                       ACHF0820
        end                                                             ACHF0821
                                                                        ACHF0822
        function gexpi(a,u,b,n)                                         ACHF0823
*                                                                       ACHF0824
*  calculates                                                           ACHF0825
*  gexpi=exp(-a) * integral from 0 to 1 {exp(a*t) * t**u * f(t)} dt ,   ACHF0826
*  where f(t) is a binomial function, expanded in terms of              ACHF0827
*  Chebyshev polinomials, f(t) = sum_{i=0}^n b_i T^{*}_i.               ACHF0828
*  See Y.L. Luke: Algorithms for the Computation of Mathematical        ACHF0829
*  Functions, Academic Press, New York, 1977, p. 126-129,               ACHF0830
*  formula (19) and (20).                                               ACHF0831
*                                                                       ACHF0832
        implicit complex*16 (a-h,o-z)                                   ACHF0833
        real*8 eps,test                                                 ACHF0834
        parameter (ldim=100,ml=1,mu=2,m=4,nband=5,eps=1.d-4)            ACHF0835
        common /iout/ iout                                              ACHF0836
        dimension b(ldim),g(ldim),band(nband,ldim),ipvt(ldim)           ACHF0837
*                                                                       ACHF0838
*  set up the system of linear band equations                           ACHF0839
*  for the Chebyshev coefficients of gexpi                              ACHF0840
*                                                                       ACHF0841
        n=min(n,ldim)                                                   ACHF0842
*                                                                       ACHF0843
   1    g(1)=2*b(1)-b(2)                                                ACHF0844
        do 2 i=2,n-1                                                    ACHF0845
          g(i)=b(i)-b(i+1)                                              ACHF0846
   2    continue                                                        ACHF0847
        g(n)=b(n)                                                       ACHF0848
        do 5 j=1,n                                                      ACHF0849
          i1=max(1,j-mu)                                                ACHF0850
          i2=min(n,j+ml)                                                ACHF0851
          do 3 i=i1,i2                                                  ACHF0852
            k=i-j+m                                                     ACHF0853
            if(i.eq.j+1) then                                           ACHF0854
              band(k,j)=a/4                                             ACHF0855
            elseif(i.eq.j) then                                         ACHF0856
              band(k,j)=a/4+u+i                                         ACHF0857
            elseif(i.eq.j-1) then                                       ACHF0858
              band(k,j)=-a/4-u+i-1                                      ACHF0859
              if(i.eq.1) band(k,j)=-u                                   ACHF0860
            elseif(i.eq.j-2) then                                       ACHF0861
              band(k,j)=-a/4                                            ACHF0862
            endif                                                       ACHF0863
   3      continue                                                      ACHF0864
   5    continue                                                        ACHF0865
*                                                                       ACHF0866
*  solution of the system linear equations                              ACHF0867
*                                                                       ACHF0868
        call zgbfa(band,nband,n,ml,mu,ipvt,info)                        ACHF0869
        if(info.ne.0) write(iout,'('' zgbfa info ='',i3)') info         ACHF0870
        call zgbsl(band,nband,n,ml,mu,ipvt,g,0)                         ACHF0871
        g(1)=g(1)/2                                                     ACHF0872
*                                                                       ACHF0873
*  gexpi=sum_{i=0}^{n-1} g(i), 1/(u+1)=sum_{i=0}^{n-1} (-1)**i*g(i)     ACHF0874
*                                                                       ACHF0875
        f0=0.d0                                                         ACHF0876
        ip=1                                                            ACHF0877
        ss=(0.d0,0.d0)                                                  ACHF0878
        do 8 i=1,n                                                      ACHF0879
          f0=f0+ip*g(i)                                                 ACHF0880
          ip=-ip                                                        ACHF0881
          ss=ss+g(i)                                                    ACHF0882
   8    continue                                                        ACHF0883
        test=abs(f0-1/(u+1))                                            ACHF0884
        if(test.gt.eps) then                                            ACHF0885
          if(n.eq.ldim) then                                            ACHF0886
            write(iout,*) ' estimated accuracy in gexpi ', test         ACHF0887
            goto 9                                                      ACHF0888
          endif                                                         ACHF0889
          n=min(n+10,ldim)                                              ACHF0890
          goto 1                                                        ACHF0891
        endif                                                           ACHF0892
   9    gexpi=ss                                                        ACHF0893
*                                                                       ACHF0894
        return                                                          ACHF0895
        end                                                             ACHF0896
                                                                        ACHF0897
        subroutine ch1f0(ap,w,nm,c)                                     ACHF0898
*                                                                       ACHF0899
*  Calculates the Chebyshev expansion of the binomial function.         ACHF0900
*  See Y.L. Luke: Algorithms for the Computation of Mathematical        ACHF0901
*  Functions, Academic Press, New York, 1977, p. 53.                    ACHF0902
*                                                                       ACHF0903
*       (1+w*t)^(-ap)=sum_{i=0}^{nm-1} c_i (w) * T^{*}_i (t)            ACHF0904
*                                                                       ACHF0905
        implicit complex*16 (a-h,o-z)                                   ACHF0906
        real*8 eps,test1                                                ACHF0907
        parameter (eps=1.d-10)                                          ACHF0908
        common /iout/ iout                                              ACHF0909
        dimension c(1)                                                  ACHF0910
*                                                                       ACHF0911
*  test wheather ap is a negative integer                               ACHF0912
*                                                                       ACHF0913
        iap=iabs(int(dble(ap)))                                         ACHF0914
        if(ap.eq.dcmplx(-iap)) then                                     ACHF0915
          n=iap-1                                                       ACHF0916
          do 1 i=n+2,nm                                                 ACHF0917
   1         c(i)=0.d0                                                  ACHF0918
        else                                                            ACHF0919
          n=nm-2                                                        ACHF0920
        endif                                                           ACHF0921
*                                                                       ACHF0922
*  computes coefficients by means of backward recurence scheme          ACHF0923
*                                                                       ACHF0924
        start=(2.d0*(1.d0+2.d0/dble(w)))**(-n)                          ACHF0925
        a2=0.d0                                                         ACHF0926
        a1=start                                                        ACHF0927
        nc=n+2                                                          ACHF0928
        c(nc)=start                                                     ACHF0929
        x1=n+2                                                          ACHF0930
        v1=1.d0-ap                                                      ACHF0931
        afac=2.d0+4.d0/w                                                ACHF0932
        do 2 k=1,n+1                                                    ACHF0933
           x1=x1-1.d0                                                   ACHF0934
           nc=nc-1                                                      ACHF0935
           c(nc)=-(x1*afac*a1+(x1+v1)*a2)/(x1-v1)                       ACHF0936
           a2=a1                                                        ACHF0937
           a1=c(nc)                                                     ACHF0938
   2    continue                                                        ACHF0939
        c(1)=c(1)/2.d0                                                  ACHF0940
        rho=c(1)                                                        ACHF0941
        p=1.d0                                                          ACHF0942
        do 3 i=2,n+2                                                    ACHF0943
           rho=rho-p*c(i)                                               ACHF0944
           p=-p                                                         ACHF0945
   3    continue                                                        ACHF0946
        do 4 i=1,n+2                                                    ACHF0947
           c(i)=c(i)/rho                                                ACHF0948
   4    continue                                                        ACHF0949
*                                                                       ACHF0950
*   accuracy test                                                       ACHF0951
*                                                                       ACHF0952
        if(n.gt.1) then                                                 ACHF0953
           sum1=-c(2)+4.d0*c(3)                                         ACHF0954
           p=-1.d0                                                      ACHF0955
           do 5 i=2,n                                                   ACHF0956
              sum1=sum1+p*(i+1)*(i+1)*c(i+2)                            ACHF0957
              p=-p                                                      ACHF0958
   5       continue                                                     ACHF0959
           test1=abs(sum1-ap*w/2.d0)                                    ACHF0960
           if(test1.gt.eps)  write(iout,'(a,d12.3)')                    ACHF0961
     &        ' estimated accuracy in ch1f0 : ',test1                   ACHF0962
        endif                                                           ACHF0963
*                                                                       ACHF0964
        return                                                          ACHF0965
        end                                                             ACHF0966
*                                                                       ACHF0967
*  further subroutines from the LINPACK library:                        ACHF0968
*     zsifa, zsidi, zsisl, zgbfa, zgbsl, zgbcal, dcabs1,                ACHF0969
*     izamax, zdotu, zaxpy, zdotc    and                                ACHF0970
*  Coulomb function RCWF by                                             ACHF0971
*   A.R.Barnett, D.H.Feng, J.W.Steed,L.J.B.Goldfarb, CPC 8 (1974) 377   ACHF0972
*--------------------------------------------------------------------   ACHF0973
*                                                                       ACHF0974
*  one may run the code in UNIX as  a.out < ao.sd > ao.sr               ACHF0975
*                                                                       ACHF0976
*  input file  ao.sd                                                    ACHF0977
*                                                                       ACHF0978
'gisp00.d','gispm05.d',                                                 ACHF0979
 4.0235d0,16.088d0,'ao16p.scd',                                         ACHF0980
 3, 1.75,                                                               ACHF0981
 24, 24,                                                                ACHF0982
 12.d0, 25, .2,                                                         ACHF0983
*                                                                       ACHF0984
*  input file gisp00.d                                                  ACHF0985
*                                                                       ACHF0986
          48                 0.                                         ACHF0987
  2.981123582996164E-002  7.426200582802977E-002                        ACHF0988
  0.157107990617875       0.152271949809351                             ACHF0989
  0.386265037576454       0.190409088263911                             ACHF0990
  0.717574694116971       0.186633059484807                             ACHF0991
   1.15139383402644       0.153424200157578                             ACHF0992
   1.68818582341905       0.108779692807490                             ACHF0993
   2.32852700665323       6.746073860921927E-002                        ACHF0994
   3.07311086165265       3.688119411582105E-002                        ACHF0995
   3.92275241304649       1.785684426915667E-002                        ACHF0996
   4.87839335592135       7.677616514497570E-003                        ACHF0997
   5.94110805462456       2.935785903739469E-003                        ACHF0998
   7.11211053589075       9.990655378158767E-004                        ACHF0999
   8.39276259909124       3.025980169922555E-004                        ACHF1000
   9.78458318468734       8.153871180355333E-005                        ACHF1001
   11.2892591680095       1.953158715728041E-005                        ACHF1002
   12.9086577782855       4.154182945052104E-006                        ACHF1003
   14.6448408832097       7.833700380277454E-007                        ACHF1004
   16.5000814289646       1.307394774920574E-007                        ACHF1005
   18.4768823868741       1.927071408017016E-008                        ACHF1006
   20.5779986340222       2.502638937126329E-009                        ACHF1007
   22.8064622905214       2.855785508771569E-010                        ACHF1008
   25.1656121564391       2.854622412059074E-011                        ACHF1009
   27.6591280444806       2.491010684937161E-012                        ACHF1010
   30.2910710010086       1.890336606971472E-013                        ACHF1011
   33.0659306624988       1.242162685949120E-014                        ACHF1012
   35.9886813274790       7.034231520212479E-016                        ACHF1013
   39.0648487641978       3.414549148591794E-017                        ACHF1014
   42.3005903629031       1.412315414895672E-018                        ACHF1015
   45.7027920385115       4.944218008097224E-020                        ACHF1016
   49.2791863828369       1.453952481367813E-021                        ACHF1017
   53.0384980878167       3.561068365003849E-023                        ACHF1018
   56.9906248148045       7.194055996494519E-025                        ACHF1019
   61.1468647861403       1.185537228350548E-026                        ACHF1020
   65.5202069290186       1.573491357075550E-028                        ACHF1021
   70.1257062361133       1.657285440919349E-030                        ACHF1022
   74.9809775189113       1.361434162716336E-032                        ACHF1023
   80.1068573503244       8.546155813962816E-035                        ACHF1024
   85.5283111160343       4.000090532480939E-037                        ACHF1025
   91.2757079936682       1.355019991102786E-039                        ACHF1026
   97.3866677135817       3.201636795354339E-042                        ACHF1027
   103.908833357176       5.035869166060774E-045                        ACHF1028
   110.904220884976       4.962487540702082E-048                        ACHF1029
   118.456425046284       2.823510716119612E-051                        ACHF1030
   126.683425768886       8.268446069503012E-055                        ACHF1031
   135.762589577865       1.049064847820770E-058                        ACHF1032
   145.986432709463       4.346574422738855E-063                        ACHF1033
   157.915612022978       3.434736438395980E-068                        ACHF1034
   172.996328148563       1.319066088398003E-074                        ACHF1035
*                                                                       ACHF1036
*  input file gispm05.d                                                 ACHF1037
*                                                                       ACHF1038
          48 -0.500000000000000                                         ACHF1039
  1.278457233711745E-002  0.446540046995086                             ACHF1040
  0.115081483588378       0.403226026958469                             ACHF1041
  0.319783878279432       0.328759261196996                             ACHF1042
  0.627109535856757       0.241966233725031                             ACHF1043
   1.03738672215406       0.160707900161284                             ACHF1044
   1.55105613079625       9.627981346217995E-002                        ACHF1045
   2.16867351479483       5.200068054248895E-002                        ACHF1046
   2.89091304669443       2.530275022560015E-002                        ACHF1047
   3.71857145710970       1.108326257887788E-002                        ACHF1048
   4.65257301434740       4.366263952193433E-003                        ACHF1049
   5.69397542244271       1.545391643363289E-003                        ACHF1050
   6.84397673183353       4.908361854610483E-004                        ACHF1051
   8.10392337664535       1.397086488908303E-004                        ACHF1052
   9.47531947589192       3.558361183644027E-005                        ACHF1053
   10.9598375637347       8.096466357554405E-006                        ACHF1054
   12.5593309474486       1.642708909566552E-006                        ACHF1055
   14.2758479323996       2.965938810009984E-007                        ACHF1056
   16.1116482030640       4.754732429514802E-008                        ACHF1057
   18.0692217104087       6.751169462255367E-009                        ACHF1058
   20.1513104920617       8.467187567022060E-010                        ACHF1059
   22.3609339469653       9.352020771970011E-011                        ACHF1060
   24.7014182063956       9.066606749167704E-012                        ACHF1061
   27.1764303961281       7.687363449753173E-013                        ACHF1062
   29.7900187807575       5.677534790338497E-014                        ACHF1063
   32.5466600353495       3.636338260986426E-015                        ACHF1064
   35.4513152220926       2.009810298315003E-016                        ACHF1065
   38.5094964891880       9.533662342373185E-018                        ACHF1066
   41.7273470970119       3.857758653674308E-019                        ACHF1067
   45.1117381723318       1.322590235296160E-020                        ACHF1068
   48.6703866831494       3.812520377375046E-022                        ACHF1069
   52.4120006467825       9.161202197099637E-024                        ACHF1070
   56.3464597342432       1.817193129409686E-025                        ACHF1071
   60.4850425305087       2.942483109948401E-027                        ACHF1072
   64.8407162572861       3.839941413589497E-029                        ACHF1073
   69.4285115890426       3.979093788538901E-031                        ACHF1074
   74.2660156887036       3.217746415876857E-033                        ACHF1075
   79.3740331847917       1.989360011844198E-035                        ACHF1076
   84.7774918932803       9.174803331715442E-038                        ACHF1077
   90.5067159166133       3.063597929422459E-040                        ACHF1078
   96.5992696619433       7.137877123560422E-043                        ACHF1079
   103.102726489260       1.107412963956401E-045                        ACHF1080
   110.079011673338       1.076644501393989E-048                        ACHF1081
   117.611597334412       6.044587863003915E-052                        ACHF1082
   125.818289122308       1.746746567183330E-055                        ACHF1083
   134.876188756134       2.186764947502585E-059                        ACHF1084
   145.077369431289       8.937404843361401E-064                        ACHF1085
   156.981623103640       6.961682775435475E-069                        ACHF1086
   172.032866745166       2.630674295402287E-075                        ACHF1087
*                                                                       ACHF1088
*  input file ao16p.scd                                                 ACHF1089
*                                                                       ACHF1090
                                                                        ACHF1091
1,16,                                                                   ACHF1092
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1093
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1094
                                                                        ACHF1095
0,                                                                      ACHF1096
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1097
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1098
   0.318154709d0,4,                                                     ACHF1099
   300.d0,                                                              ACHF1100
   300.d0,                                                              ACHF1101
   300.d0,                                                              ACHF1102
   300.d0,                                                              ACHF1103
   0.7708d0,                                                            ACHF1104
   0.4897d0,                                                            ACHF1105
   0.2815d0,                                                            ACHF1106
   0.1541d0,                                                            ACHF1107
   0.0822d0,                                                            ACHF1108
   0.0432d0,                                                            ACHF1109
   0.0225d0,                                                            ACHF1110
   0.0116d0,                                                            ACHF1111
   0.0059d0,                                                            ACHF1112
   0.0030d0,                                                            ACHF1113
   0.0015d0,                                                            ACHF1114
   0.0008d0,                                                            ACHF1115
   0.0004d0,                                                            ACHF1116
   0.0002d0,                                                            ACHF1117
   0.0001d0,                                                            ACHF1118
   0.d0,                                                                ACHF1119
                                                                        ACHF1120
1,                                                                      ACHF1121
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1122
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1123
   0.318154709d0,4,                                                     ACHF1124
   300.d0,                                                              ACHF1125
   300.d0,                                                              ACHF1126
   300.d0,                                                              ACHF1127
   300.d0,                                                              ACHF1128
   0.6562d0,                                                            ACHF1129
   0.3804d0,                                                            ACHF1130
   0.21d0,                                                              ACHF1131
   0.1129d0,                                                            ACHF1132
   0.0597d0,                                                            ACHF1133
   0.0312d0,                                                            ACHF1134
   0.0161d0,                                                            ACHF1135
   0.0083d0,                                                            ACHF1136
   0.0042d0,                                                            ACHF1137
   0.0021d0,                                                            ACHF1138
   0.0011d0,                                                            ACHF1139
   0.0005d0,                                                            ACHF1140
   0.0002d0,                                                            ACHF1141
   0.0001d0,                                                            ACHF1142
   0.d0,                                                                ACHF1143
                                                                        ACHF1144
2,                                                                      ACHF1145
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1146
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1147
   0.318154709d0,3,                                                     ACHF1148
   300.d0,                                                              ACHF1149
   300.d0,                                                              ACHF1150
   300.d0,                                                              ACHF1151
   0.7708d0,                                                            ACHF1152
   0.4897d0,                                                            ACHF1153
   0.2815d0,                                                            ACHF1154
   0.1541d0,                                                            ACHF1155
   0.0822d0,                                                            ACHF1156
   0.0432d0,                                                            ACHF1157
   0.0225d0,                                                            ACHF1158
   0.0116d0,                                                            ACHF1159
   0.0059d0,                                                            ACHF1160
   0.0030d0,                                                            ACHF1161
   0.0015d0,                                                            ACHF1162
   0.0008d0,                                                            ACHF1163
   0.0004d0,                                                            ACHF1164
   0.0002d0,                                                            ACHF1165
   0.0001d0,                                                            ACHF1166
   0.d0,                                                                ACHF1167
                                                                        ACHF1168
3,                                                                      ACHF1169
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1170
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1171
   0.318154709d0,3,                                                     ACHF1172
   300.d0,                                                              ACHF1173
   300.d0,                                                              ACHF1174
   300.d0,                                                              ACHF1175
   0.6562d0,                                                            ACHF1176
   0.3804d0,                                                            ACHF1177
   0.21d0,                                                              ACHF1178
   0.1129d0,                                                            ACHF1179
   0.0597d0,                                                            ACHF1180
   0.0312d0,                                                            ACHF1181
   0.0161d0,                                                            ACHF1182
   0.0083d0,                                                            ACHF1183
   0.0042d0,                                                            ACHF1184
   0.0021d0,                                                            ACHF1185
   0.0011d0,                                                            ACHF1186
   0.0005d0,                                                            ACHF1187
   0.0002d0,                                                            ACHF1188
   0.0001d0,                                                            ACHF1189
   0.d0,                                                                ACHF1190
                                                                        ACHF1191
4,                                                                      ACHF1192
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1193
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1194
   0.318154709d0,2,                                                     ACHF1195
   300.d0,                                                              ACHF1196
   300.d0,                                                              ACHF1197
   0.7708d0,                                                            ACHF1198
   0.4897d0,                                                            ACHF1199
   0.2815d0,                                                            ACHF1200
   0.1541d0,                                                            ACHF1201
   0.0822d0,                                                            ACHF1202
   0.0432d0,                                                            ACHF1203
   0.0225d0,                                                            ACHF1204
   0.0116d0,                                                            ACHF1205
   0.0059d0,                                                            ACHF1206
   0.0030d0,                                                            ACHF1207
   0.0015d0,                                                            ACHF1208
   0.0008d0,                                                            ACHF1209
   0.0004d0,                                                            ACHF1210
   0.0002d0,                                                            ACHF1211
   0.0001d0,                                                            ACHF1212
   0.d0,                                                                ACHF1213
                                                                        ACHF1214
5,                                                                      ACHF1215
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1216
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1217
   0.318154709d0,2,                                                     ACHF1218
   300.d0,                                                              ACHF1219
   300.d0,                                                              ACHF1220
   0.6562d0,                                                            ACHF1221
   0.3804d0,                                                            ACHF1222
   0.21d0,                                                              ACHF1223
   0.1129d0,                                                            ACHF1224
   0.0597d0,                                                            ACHF1225
   0.0312d0,                                                            ACHF1226
   0.0161d0,                                                            ACHF1227
   0.0083d0,                                                            ACHF1228
   0.0042d0,                                                            ACHF1229
   0.0021d0,                                                            ACHF1230
   0.0011d0,                                                            ACHF1231
   0.0005d0,                                                            ACHF1232
   0.0002d0,                                                            ACHF1233
   0.0001d0,                                                            ACHF1234
   0.d0,                                                                ACHF1235
                                                                        ACHF1236
6,                                                                      ACHF1237
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1238
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1239
   0.318154709d0,1,                                                     ACHF1240
   300.d0,                                                              ACHF1241
   0.7708d0,                                                            ACHF1242
   0.4897d0,                                                            ACHF1243
   0.2815d0,                                                            ACHF1244
   0.1541d0,                                                            ACHF1245
   0.0822d0,                                                            ACHF1246
   0.0432d0,                                                            ACHF1247
   0.0225d0,                                                            ACHF1248
   0.0116d0,                                                            ACHF1249
   0.0059d0,                                                            ACHF1250
   0.0030d0,                                                            ACHF1251
   0.0015d0,                                                            ACHF1252
   0.0008d0,                                                            ACHF1253
   0.0004d0,                                                            ACHF1254
   0.0002d0,                                                            ACHF1255
   0.0001d0,                                                            ACHF1256
   0.d0,                                                                ACHF1257
                                                                        ACHF1258
7,                                                                      ACHF1259
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1260
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1261
   0.318154709d0,1,                                                     ACHF1262
   300.d0,                                                              ACHF1263
   0.6562d0,                                                            ACHF1264
   0.3804d0,                                                            ACHF1265
   0.21d0,                                                              ACHF1266
   0.1129d0,                                                            ACHF1267
   0.0597d0,                                                            ACHF1268
   0.0312d0,                                                            ACHF1269
   0.0161d0,                                                            ACHF1270
   0.0083d0,                                                            ACHF1271
   0.0042d0,                                                            ACHF1272
   0.0021d0,                                                            ACHF1273
   0.0011d0,                                                            ACHF1274
   0.0005d0,                                                            ACHF1275
   0.0002d0,                                                            ACHF1276
   0.0001d0,                                                            ACHF1277
   0.d0,                                                                ACHF1278
                                                                        ACHF1279
8,                                                                      ACHF1280
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1281
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1282
   0.318154709d0,0,                                                     ACHF1283
   0.7708d0,                                                            ACHF1284
   0.4897d0,                                                            ACHF1285
   0.2815d0,                                                            ACHF1286
   0.1541d0,                                                            ACHF1287
   0.0822d0,                                                            ACHF1288
   0.0432d0,                                                            ACHF1289
   0.0225d0,                                                            ACHF1290
   0.0116d0,                                                            ACHF1291
   0.0059d0,                                                            ACHF1292
   0.0030d0,                                                            ACHF1293
   0.0015d0,                                                            ACHF1294
   0.0008d0,                                                            ACHF1295
   0.0004d0,                                                            ACHF1296
   0.0002d0,                                                            ACHF1297
   0.0001d0,                                                            ACHF1298
   0.d0,                                                                ACHF1299
                                                                        ACHF1300
9,                                                                      ACHF1301
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1302
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1303
   0.318154709d0,0,                                                     ACHF1304
   0.6562d0,                                                            ACHF1305
   0.3804d0,                                                            ACHF1306
   0.21d0,                                                              ACHF1307
   0.1129d0,                                                            ACHF1308
   0.0597d0,                                                            ACHF1309
   0.0312d0,                                                            ACHF1310
   0.0161d0,                                                            ACHF1311
   0.0083d0,                                                            ACHF1312
   0.0042d0,                                                            ACHF1313
   0.0021d0,                                                            ACHF1314
   0.0011d0,                                                            ACHF1315
   0.0005d0,                                                            ACHF1316
   0.0002d0,                                                            ACHF1317
   0.0001d0,                                                            ACHF1318
   0.d0,                                                                ACHF1319
                                                                        ACHF1320
10,                                                                     ACHF1321
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1322
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1323
   0.318154709d0,0,                                                     ACHF1324
   0.4897d0,                                                            ACHF1325
   0.2815d0,                                                            ACHF1326
   0.1541d0,                                                            ACHF1327
   0.0822d0,                                                            ACHF1328
   0.0432d0,                                                            ACHF1329
   0.0225d0,                                                            ACHF1330
   0.0116d0,                                                            ACHF1331
   0.0059d0,                                                            ACHF1332
   0.0030d0,                                                            ACHF1333
   0.0015d0,                                                            ACHF1334
   0.0008d0,                                                            ACHF1335
   0.0004d0,                                                            ACHF1336
   0.0002d0,                                                            ACHF1337
   0.0001d0,                                                            ACHF1338
   0.d0,                                                                ACHF1339
                                                                        ACHF1340
11,                                                                     ACHF1341
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1342
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1343
   0.318154709d0,0,                                                     ACHF1344
   0.3804d0,                                                            ACHF1345
   0.21d0,                                                              ACHF1346
   0.1129d0,                                                            ACHF1347
   0.0597d0,                                                            ACHF1348
   0.0312d0,                                                            ACHF1349
   0.0161d0,                                                            ACHF1350
   0.0083d0,                                                            ACHF1351
   0.0042d0,                                                            ACHF1352
   0.0021d0,                                                            ACHF1353
   0.0011d0,                                                            ACHF1354
   0.0005d0,                                                            ACHF1355
   0.0002d0,                                                            ACHF1356
   0.0001d0,                                                            ACHF1357
   0.d0,                                                                ACHF1358
                                                                        ACHF1359
12,                                                                     ACHF1360
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1361
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1362
   0.318154709d0,0,                                                     ACHF1363
   0.2815d0,                                                            ACHF1364
   0.1541d0,                                                            ACHF1365
   0.0822d0,                                                            ACHF1366
   0.0432d0,                                                            ACHF1367
   0.0225d0,                                                            ACHF1368
   0.0116d0,                                                            ACHF1369
   0.0059d0,                                                            ACHF1370
   0.0030d0,                                                            ACHF1371
   0.0015d0,                                                            ACHF1372
   0.0008d0,                                                            ACHF1373
   0.0004d0,                                                            ACHF1374
   0.0002d0,                                                            ACHF1375
   0.0001d0,                                                            ACHF1376
   0.d0,                                                                ACHF1377
                                                                        ACHF1378
13,                                                                     ACHF1379
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1380
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1381
   0.318154709d0,0,                                                     ACHF1382
   0.21d0,                                                              ACHF1383
   0.1129d0,                                                            ACHF1384
   0.0597d0,                                                            ACHF1385
   0.0312d0,                                                            ACHF1386
   0.0161d0,                                                            ACHF1387
   0.0083d0,                                                            ACHF1388
   0.0042d0,                                                            ACHF1389
   0.0021d0,                                                            ACHF1390
   0.0011d0,                                                            ACHF1391
   0.0005d0,                                                            ACHF1392
   0.0002d0,                                                            ACHF1393
   0.0001d0,                                                            ACHF1394
   0.d0,                                                                ACHF1395
                                                                        ACHF1396
14,                                                                     ACHF1397
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1398
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1399
   0.318154709d0,0,                                                     ACHF1400
   0.1541d0,                                                            ACHF1401
   0.0822d0,                                                            ACHF1402
   0.0432d0,                                                            ACHF1403
   0.0225d0,                                                            ACHF1404
   0.0116d0,                                                            ACHF1405
   0.0059d0,                                                            ACHF1406
   0.0030d0,                                                            ACHF1407
   0.0015d0,                                                            ACHF1408
   0.0008d0,                                                            ACHF1409
   0.0004d0,                                                            ACHF1410
   0.0002d0,                                                            ACHF1411
   0.0001d0,                                                            ACHF1412
   0.d0,                                                                ACHF1413
                                                                        ACHF1414
15,                                                                     ACHF1415
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1416
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1417
   0.318154709d0,0,                                                     ACHF1418
   0.1129d0,                                                            ACHF1419
   0.0597d0,                                                            ACHF1420
   0.0312d0,                                                            ACHF1421
   0.0161d0,                                                            ACHF1422
   0.0083d0,                                                            ACHF1423
   0.0042d0,                                                            ACHF1424
   0.0021d0,                                                            ACHF1425
   0.0011d0,                                                            ACHF1426
   0.0005d0,                                                            ACHF1427
   0.0002d0,                                                            ACHF1428
   0.0001d0,                                                            ACHF1429
   0.d0,                                                                ACHF1430
                                                                        ACHF1431
16,                                                                     ACHF1432
 16.d0,  -139.16d0, 0.04111d0, 0.13428d0,                               ACHF1433
      -1.67709d0,  0.18963d0, 0.43546d0,                                ACHF1434
   0.318154709d0,0,                                                     ACHF1435
   0.0822d0,                                                            ACHF1436
   0.0432d0,                                                            ACHF1437
   0.0225d0,                                                            ACHF1438
   0.0116d0,                                                            ACHF1439
   0.0059d0,                                                            ACHF1440
   0.0030d0,                                                            ACHF1441
   0.0015d0,                                                            ACHF1442
   0.0008d0,                                                            ACHF1443
   0.0004d0,                                                            ACHF1444
   0.0002d0,                                                            ACHF1445
   0.0001d0,                                                            ACHF1446
   0.d0,                                                                ACHF1447
                                                                        ACHF1448
*                                                                       ACHF1449
*  output file ao.sr                                                    ACHF1450
*                                                                       ACHF1451
                                                                        ACHF1452
  Scattering States in Coulomb-like potential                           ACHF1453
                                                                        ACHF1454
                                                                        ACHF1455
  mass_1 =  4.02350000 mass_2 = 16.08800000                             ACHF1456
                                                                        ACHF1457
  Potential file name : ao16p.scd                                       ACHF1458
                                                                        ACHF1459
  lm =  3   omcst = 1.7500   b = 3.1396                                 ACHF1460
                                                                        ACHF1461
                                                                        ACHF1462
  Initial energy =      12.00                                           ACHF1463
                                                                        ACHF1464
                                                                        ACHF1465
  nmax=   24                                                            ACHF1466
                                                                        ACHF1467
  Coulomb phase shift eta_l   =    1.667727 rad                         ACHF1468
                      eta_l   =   95.553724 deg                         ACHF1469
                                                                        ACHF1470
  Colomb modified nuclear                                               ACHF1471
          phase shift delta_l =    1.508089 rad                         ACHF1472
                      delta_l =   86.407137 deg                         ACHF1473
                                                                        ACHF1474
  Coulomb modified nuclear                                              ACHF1475
  scattering amplitude   a'_l =    0.960157D-01 +i*  -0.727880D+00      ACHF1476
                                                                        ACHF1477
  wave function                                                         ACHF1478
                                                                        ACHF1479
    0.20 -0.2045D-03 -0.3257D-02                                        ACHF1480
    0.40 -0.2825D-02 -0.4500D-01                                        ACHF1481
    0.60 -0.1113D-01 -0.1773D+00                                        ACHF1482
    0.80 -0.2436D-01 -0.3879D+00                                        ACHF1483
    1.00 -0.3560D-01 -0.5670D+00                                        ACHF1484
    1.20 -0.3586D-01 -0.5711D+00                                        ACHF1485
    1.40 -0.2112D-01 -0.3364D+00                                        ACHF1486
    1.60  0.3525D-02  0.5615D-01                                        ACHF1487
    1.80  0.2618D-01  0.4170D+00                                        ACHF1488
    2.00  0.3527D-01  0.5617D+00                                        ACHF1489
    2.20  0.2631D-01  0.4190D+00                                        ACHF1490
    2.40  0.4077D-02  0.6493D-01                                        ACHF1491
    2.60 -0.2060D-01 -0.3282D+00                                        ACHF1492
    2.80 -0.3675D-01 -0.5853D+00                                        ACHF1493
    3.00 -0.3805D-01 -0.6060D+00                                        ACHF1494
    3.20 -0.2459D-01 -0.3916D+00                                        ACHF1495
    3.40 -0.1403D-02 -0.2234D-01                                        ACHF1496
    3.60  0.2445D-01  0.3894D+00                                        ACHF1497
    3.80  0.4671D-01  0.7439D+00                                        ACHF1498
    4.00  0.6146D-01  0.9789D+00                                        ACHF1499
    4.20  0.6743D-01  0.1074D+01                                        ACHF1500
    4.40  0.6530D-01  0.1040D+01                                        ACHF1501
    4.60  0.5683D-01  0.9051D+00                                        ACHF1502
    4.80  0.4402D-01  0.7011D+00                                        ACHF1503
    5.00  0.2869D-01  0.4569D+00                                        ACHF1504
                                                                        ACHF1505
                                                                        ACHF****
