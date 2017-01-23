      subroutine massdz33 (Z,X,Y)  ! Z protons, X neutrons , Y Bind.Energy
      common/cpair/epair,idef
      DIMENSIONA(33),DYDA(33),FYDA(33),op(2,3,2),fyd0(33),onps(2),n4(2),
     1ONP(0:8,2,2),OT(0:8,2,2),OEI(2),DEI(2),NN(2),NOC(18,2),OP2(2),
     2JUP(2),JUD(2),YM(2),OP1(2),n2(2),shell(2),sshell(2)
CC*****************************************
      data a/ 9.0914,   6.3355,   4.5791,  19.8946,   1.7325,   7.5247,
     &       -7.1953, -39.9787,  -0.3976,   0.8131,  -0.7435,  -3.7291,
     &       -0.1305,  -0.6387,   0.4534,   2.0605,   0.3449,   1.4727,
     &       -1.0433,   0.0000,   5.2495,   0.0000, -32.1007,-151.1164,
     &       -4.6103, -32.4238, -37.3226, -52.1673,   0.9597,   3.0024,
     &        0.6977,   6.0390,  17.7960/
      double precision Y
      integer Z,X
c Data=1751 RMS= 0.330  (mars 95)
cFM+*   9.09 fm+*   6.34 FS+*   4.58 fs+*  19.89 FS-*   1.73 fs-*   7.52
cFC+*  -7.20 fc+* -39.98 PM+*  -0.40 pm+*   0.81 PS+*  -0.74 ps+*  -3.73
cPS-*  -0.13 ps-*  -0.64 S3 *   0.45 s3 *   2.06 SQ-*   0.34 sq-*   1.47
cD3 *  -1.04 d3     0.00 QQ+*   5.25 qq+    0.00  D0* -32.10  d0*-151.12
cQQ-*  -4.61 qq-* -32.42  TT* -37.32  tt* -52.17  SS*   0.96  ss*   3.00
c C *   0.70 P0 *   6.04 P1 *  17.80
c                    ----------------------
CC DEBUT DU SS PROG.
      idef=0
      IMAX=18
      MO=2
      MAXP=8
      NN(1)=X
      NN(2)=Z
      nx=x
      nz=z
      V=X+Z
      T=ABS(X-Z)
      r=v**(1./3.)
      s=r*r
        Rc=R*(1.-.25*(t/v)**2)
        Ra=(rc*rc)/r
C******************
      DO NDEF=1,2               !  NDEF=1 sph.  NDEF=2  def.
       ym(ndef)=0.
C   NDEF= 1->SPH// 2->DEF//
      JU=0
      JUP(1)=0
      JUP(2)=0
      JUD(1)=0
      JUD(2)=0
C-----------
      IF(NDEF.EQ.2.and.nz.gt.50)JU=4
      DO KK=1,33
         FYDA(KK)=0.
         FYD0(KK)=0.
         dyda(KK)=0.
      END DO
      DO J=1,2
        DO I=1,IMAX
        NOC(I,J)=0
        ENDDO
        DO K=0,MAXP
       DO M=1,MO
        ONP(K,M,J)=0.
        ENDDO
       ENDDO
      ENDDO
      DO J=1,2                  !LOOP OVER FLUIDS (1=N,2=P)
      NCUM=0
      I=0
  20  I=I+1
        i2=(I/2)*2
         if(i2.ne.i)then
         ID=I+1
       if(ncum.lt.nn(j))sshell(J)=1.     !sscouche  J
         else
         ID=I*(I-2)/4
       if(ncum.lt.nn(j))sshell(J)=2.     ! SSC R
         endif
      NCUM=NCUM+ID
      IF(NCUM.LE.NN(J))THEN
      NOC(I,J)=ID
      GO TO 20
      ENDIF
          shell(J)=I        !N0 de sscouche  (SSc 2-> 0 nucl)
          IP=(I-1)/2        !N0 couche HO
      MOC=NN(J)-NCUM+ID
      IF(NDEF.EQ.2) THEN
         if(i2.ne.i)then
      JUD(J)=MAX (0,JU-MOC)
      JUP(J)=0
      Else
      JUP(J)=MIN (JU,MOC)
      JUP(J)=JU
      JUD(J)=0
      ENDIF
      ENDIF
      NOC(I,J)=MOC-JUP(J)+JUD(J)
      NOC(I+1,J)=JUP(j)
      NOC(I-1,J)=NOC(I-1,J)-JUD(J)
         if(i2.ne.i)then
      OEI(J)=MOC+IP*(IP-1)-JU
      DEI(J)=IP*(IP+1)+2
      ELSE
      OEI(J)=MOC-JU
      DEI(J)=(IP+1)*(IP+2)+2
      ENDIF
CC HERE,DEGENERACIES AND NUMBER OF ACTIVE PARTICLES  FOR  EI.
      IPL=0
         vmr=0.
         vmj=0.
      DO II=1,IMAX
         onps(j)=0.
         IP=(II-1)/2
         DEGI=(IP+1)*(IP+2)
         FAC=1./SQRT(DEGI)
      IF(IP.NE.IPL)IPL=IPL+1
        if((2*IP+1).eq.II)then
        VM2=(.5*IP)/(IP+1)
        degr=ip*(ip-1)
         if(ip.gt.2)then
        vmr=(.5*(ip-1))/ip
        vmj=-1./ip
         if(noc(ii,j).le.degr)onps(j)=noc(ii,j)*vmr
         if(noc(ii,j).gt.degr)onps(j)=degr*vmr+(noc(ii,j)-degr)*vmj
            endif
        endif
        if((2*IP+1).ne.II)VM2=-1./(IP+1)      !SSc. j
      ONP(IPL,2,J)=ONP(IPL,2,J)+NOC(II,J)*VM2
      ONP(IPL,1,J)=ONP(IPL,1,J)+NOC(II,J)*FAC
            fyd0(29)=fyd0(29)+onps(j)*(onp(ipl,1,j)+onp(ipl,2,j))
      ENDDO
      ENDDO           !END OF LOOP OVER FLUIDS
c**************
      IF(NDEF.EQ.2) THEN
      ALFA=0.
      ELSE
      ALFA=0.5
      ENDIF
      FACN=DEI(1)**ALFA
      FACZ=DEI(2)**ALFA
      DNNB=OEI(1)*(DEI(1)-OEI(1))/DEI(1)
      DZZB=OEI(2)*(DEI(2)-OEI(2))/DEI(2)
      QN=DNNB*FACN/SQRT(DEI(1))
      QZ=DZZB*FACZ/SQRT(DEI(2))
      DI1N=DNNB*(2*OEI(1)-DEI(1))*FACN*FACN/(DEI(1))
      DI1Z=DZZB*(2*OEI(2)-DEI(2))*FACZ*FACZ/(DEI(2))
      S3=DI1Z+DI1N
      QQ0=(QN+QZ)*(QN+QZ)
      QQ1=(QN-QZ)*(QN-QZ)
      qqp=qq0+qq1
      qqm=qq0-qq1
c------------
      DO M=1,MO
      DO I=0,MAXP
      OT(I,M,1)=ONP(I,M,1)+ONP(I,M,2)
      OT(I,M,2)=ONP(I,M,1)-ONP(I,M,2)
      ENDDO
      ENDDO
c------------
      DO L=1,2
      DO M=1,3
      DO J=1,2
      OP(L,M,J)=0.
      ENDDO
      ENDDO
      ENDDO
      DO I=0,MAXP
      DEGI=(I+1)*(I+2)
      FAC=SQRT(DEGI)
      DO N=1,2
      DO M=1,MO
      OP(1,M,N)=OP(1,M,N)+OT(I,M,N)
       OTX=OT(I,M,N)*OT(I,M,N)*FAC
      OP(2,M,N)=OP(2,M,N)+OTX
      ENDDO
      ENDDO
      ENDDO
      DO N=1,2
       OP1(N)=0.
       OP2(N)=0.
      DO M=1,MO
      OPXX=OP(1,M,N)*OP(1,M,N)
      OP(1,M,N)=OPXX
      ENDDO
      ENDDO
      DO N=1,2
      DO I=0,MAXP
      DEGI=(I+1)*(I+2)
      FAC=sqrt(degi)
      OP1(N)=OP1(N)+OT(I,1,N)/FAC
      OP2(N)=OP2(N)+OT(I,2,N)/FAC
      ENDDO
      ENDDO
      DO N=1,2
      OP(1,3,N)=OP1(N)*OP2(N)
      ENDDO
      K=-1
      DO L=1,2
      DO M=1,3
      DO N=1,2
      K=K+2
      FYD0(K)=OP(L,M,N)
      ENDDO
      ENDDO
      ENDDO
C-----------
        DO JF=1,17,4
        FYD0(JF)=FYD0(JF)+FYD0(JF+2)
        FYD0(JF+2)=FYD0(JF)-2.*FYD0(JF+2)
        ENDDO
C-----------
        fyda(1)=fyd0(1)
        fyda(3)=fyd0(5)
        fyda(5)=fyd0(7)
        fyda(7)=fyd0(9)
        fyda(9)=fyd0(13)
        fyda(11)=fyd0(17)
        fyda(13)=fyd0(19)
      FYDA(27)= T*(T+2)/s
      fyda(29)=fyd0(29)
C-----------
      IF(NDEF.EQ.1)then
       FYDA(15)=S3
      FYDA(17)=QQm
      ELSE
      FYDA(19)=S3
      FYDA(21)=QQp
      FYDA(23)=16.-qqm
      FYDA(25)=QQm
      ENDIF
C-----------
      DO MSS=1,29,2
       DYDA(MSS)=FYDA(MSS)/Ra
       DYDA(MSS+1)=-DYDA(MSS)/Ra
       ym(ndef)=ym(ndef)+dyda(mss)*a(mss)+dyda(mss+1)*a(mss+1)
      ENDDO
c---------
       Z2=Z*(Z-1)
      DYDA(31)=(-Z2+.76*Z2**(2./3.))/Rc      !Coulomb
c---------
         rxz=1./Ra                           !Pairing
         vxz=1./v
         txz=rxz*(t/v)
         uxz=rxz-txz
       do l=1,2
      n2(l)=2*(nn(l)/2)
      dyda(32)=-rxz+dyda(32)
      if(sshell(l).eq.2.) dyda(33)=vxz+dyda(33)
       !effet de couche en 1/A
      if(n2(l).eq.nn(l))dyda(32)=uxz+dyda(32)
         enddo
      j=2                     !  Z>N
      if(nn(1).ge.nn(2))j=1   !  N>ou=Z
      k=3-j
      if(n2(j).eq.nn(j).and.n2(k).ne.nn(k))dyda(32)=-txz+dyda(32)
c********
      DO MSS=31,33
       ym(ndef)=ym(ndef)+dyda(mss)*a(mss)
      ENDDO
      ENDDO  !End sph and def calculations
        Y=ym(1)
        if(z.gt.50.and.ym(2).gt.ym(1)) Y=ym(2)
        if(z.gt.50.and.ym(2).gt.ym(1)) idef=1
cc	av=1.717
cc 	v4=v/(1-(t/2./v)**2.)**2*(1.+2./9.*(t/v)**2)
c       v4=v
c       Y=Y-v4*(a(1)+a(9))*av+v4/Ra*(a(2)+a(10))*av
c       Y=Y-dyda(1)*a(1)+dyda(2)*a(2)-dyda(9)*a(9)+dyda(10)*a(10)
cc        t22=t*(t+2)
cc        Y=Y-15.07*v**(4./3.)/Ra+15.00*v**(4./3.)/Ra**2
cc     &    +33.88*t22/v**(2./3.)/Ra-48.51*t22/v**(2./3.)/Ra**2
        epair=dyda(32)*a(32)+dyda(33)*a(33)
      END

      subroutine massdz10(nz,nx,E)     ! Duflo-Zuker fevrier 1996                 
c Calculation of binding energy E (nx neutrons,nz protons)                      
      dimension b(10),dyda(10),op(2),n2(2),dx(2),qx(2),os(2),                   
     &          onp(0:8,2,2),oei(2),dei(2),nn(2),noc(18,2),pp(2),y(2)           
      data b/0.7043,17.7418,16.2562,37.5562,53.9017,0.4711,2.1307,              
     &       0.0210,40.5356,6.0632/
      double precision E
      integer nz,nx                                               
c*********                                                                      
      nn(1)=nx                                                                  
      nn(2)=nz                                                                  
      a=nx+nz                                                                   
      t=abs(nx-nz)                                                              
      r=a**(1./3.)                                                              
      s=r*r                                                                     
      rc=r*(1.-.25*(t/a)**2)       !      Charge radius                         
      ra=(rc*rc)/r                                                              
c--------                                                                       
      z2=nz*(nz-1)                                                              
      dyda(1)=(-z2+.76*z2**(2./3.))/rc  ! Coulomb energy                        
c********                          ! beginning of main loop                     
      do ndef=1,2                  !      ndef=1  spherical                     
      ju=0                         !      ndef=2  deformed                      
      y(ndef)=0.                                                                
      if(ndef.eq.2) ju=4           !      nucleons associated to deform.        
      do kk=2,10                                                                
        dyda(kk)=0.                                                             
      enddo                                                                     
c--------                          ! beginning of loop over N and Z             
      do j=1,2                                                                  
        do l=1,18                                                               
          noc(l,j)=0                                                            
        enddo                                                                   
        do l=1,2                                                                
          do k=0,8                                                              
            onp(k,l,j)=0.                                                       
          enddo                                                                 
        enddo                                                                   
        n2(j)=2*(nn(j)/2)          !      (for pairing calculation)             
        ncum=0                                                                  
        i=0                                                                     
c--------                                                                       
  20    i=i+1                      !     sub-shells (ssh) j and r filling       
        i2=(i/2)*2                                                              
        if(i2.ne.i)then                                                         
          id=i+1                   !             for ssh j                      
        else                                                                    
          id=i*(i-2)/4             !             for ssc r                      
        endif                                                                   
        ncum=ncum+id                                                            
        if(ncum.lt.nn(j))then                                                   
          noc(i,j)=id              !     nb of nucleons in each ssh             
          go to 20                                                              
        endif                                                                   
c--------                                                                       
        imax=i+1                   !     imax = last subshell nb                
        ip=(i-1)/2                 !     HO number (p)                          
        ipm=i/2                                                                 
        pp(j)=ip                                                                
        moc=nn(j)-ncum+id                                                       
        noc(i,j)=moc-ju            !     nb of nucleons in last ssh             
        noc(i+1,j)=ju                                                           
        if(i2.ne.i)then            !     ssh j                                  
          oei(j)=moc+ip*(ip-1)     !       nb of nucleons in last EI shell      
          dei(j)=ip*(ip+1)+2       !       size of the EI shell                 
        else                       !     ssh r                                  
          oei(j)=moc-ju            !       nb of nucleons in last EI shell      
          dei(j)=(ip+1)*(ip+2)+2   !       size of the EI shell                 
        endif                                                                   
        qx(j)=oei(j)*(dei(j)-oei(j)-ju)/dei(j)  ! n*(D-n)/D        S3(j)        
        dx(j)=qx(j)*(2*oei(j)-dei(j))           ! n*(D-n)*(2n-D)/D  Q           
        if(ndef.eq.2)qx(j)=qx(j)/sqrt(dei(j))   ! scaling for deformed          
c--------                                                                       
        do i=1,imax                             ! Amplitudes                    
          ip=(i-1)/2                                                            
          fact=sqrt((ip+1.)*(ip+2.))                                            
          onp(ip,1,j)=onp(ip,1,j)+noc(i,j)/fact !    for FM term                
          vm=-1.                                                                
          if((2*(i/2)).ne.i)vm=.5*ip            !    for spin-orbit term        
          onp(ip,2,j)=onp(ip,2,j)+noc(i,j)*vm                                   
        enddo                                                                   
c--------                                                                       
        op(j)=0.                                                                
        os(j)=0.                                                                
        do ip=0,ipm                !       FM and SO terms                      
          pi=ip                                                                 
          den=((pi+1)*(pi+2))**(3./2.)                                          
          op(j)=op(j)+onp(ip,1,j)                                ! FM           
          os(j)=os(j)+onp(ip,2,j)*(1.+onp(ip,1,j))*(pi*pi/den)   ! SO           
     &               +onp(ip,2,j)*(1.-onp(ip,1,j))*((4*pi-5)/den)               
        enddo                                                                   
        op(j)=op(j)*op(j)                                                       
      enddo                                                                     
c--------                          ! end of loop over  N and Z                  
      dyda(2)=op(1)+op(2)                 !   Master term (FM): volume          
      dyda(3)=-dyda(2)/ra                 !                     surface         
      dyda(2)=dyda(2)+os(1)+os(2)         !   FM + SO                           
      dyda(4)=-t*(t+2)/(r*r)              !   isospin term : volume             
      dyda(5)=-dyda(4)/ra                 !                : surface            
      if(ndef.eq.1)then                   ! sph.                                
        dyda(6)=dx(1)+dx(2)               !   S3  volume                        
        dyda(7)=-dyda(6)/ra               !       surface                       
        px=sqrt(pp(1))+sqrt(pp(2))                                              
        dyda(8)=qx(1)*qx(2)*(2**px)       !   QQ sph.                           
      else                                ! def.                                
        dyda(9)=qx(1)*qx(2)               !   QQ deform.                        
      endif                                                                     
      dyda(5)=t*(1-t)/(a*ra**3)+dyda(5)   !   "Wigner term"                     
c--------                                 !   PAIRING                           
      if(n2(1).ne.nn(1).and.n2(2).ne.nn(2))dyda(10)= t/a                        
      if(nx.gt.nz)then                                                          
        if(n2(1).eq.nn(1).and.n2(2).ne.nn(2))dyda(10)= 1-t/a                    
        if(n2(1).ne.nn(1).and.n2(2).eq.nn(2))dyda(10)= 1                        
      else                                                                      
        if(n2(1).eq.nn(1).and.n2(2).ne.nn(2))dyda(10)= 1                        
        if(n2(1).ne.nn(1).and.n2(2).eq.nn(2))dyda(10)= 1-t/a                    
      endif                                                                     
      if(n2(2).eq.nn(2).and.n2(1).eq.nn(1))dyda(10)= 2-t/a                      
c--------                                                                       
      do mss=2,10                                                               
        dyda(mss)=dyda(mss)/ra                                                  
      enddo                                                                     
      do mss=1,10                                                               
        y(ndef)=y(ndef)+dyda(mss)*b(mss)                                        
      enddo                                                                     
c--------                            ! end of main loop                         
      enddo                                                                     
      de=y(2)-y(1)                                                              
      E=y(2)                         ! Binding Energy for def. nuclides         
      if(de.le.0..or.nz.le.50)E=y(1) !                spherical nuclides        
      return                                                                    
      end
