**********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                           MODEL                                     *
C                                                                     *
C    This file contains Fortran subroutines into which the user       *
C    must enter the relevant model equations and constants.           *
C    Consult the User's Guide for details concerning the format for   *
C    entered equations and definition of symbols.                     *
C                                                                     *
C       1. Symbol-  Parameter symbols and model constants             *
C       2. DiffEq-  System differential equations                     *
C       3. Output-  System output equations                           *
C       4. Varmod-  Error variance model equations                    *
C       5. Covmod-  Covariate model equations (ITS,MLEM)              *
C       6. Popinit- Population parameter initial values (ITS,MLEM)    *
C       7. Prior -  Parameter mean and covariance values (ID,NPD,STS) *
C       8. Sparam-  Secondary parameters                              *
C       9. Amat  -  System state matrix                               *
C                                                                     *
C**********************************************************************

C######################################################################C

        Subroutine SYMBOL
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

CC
C----------------------------------------------------------------------C
C   Enter as Indicated                                                 C
C----c-----------------------------------------------------------------C

      NDEqs   =  23   ! Enter # of Diff. Eqs.
      NSParam =  17   ! Enter # of System Parameters.
      NVparam =  2  ! Enter # of Variance Parameters.
      NSecPar =  0  ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' interaction '

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C

	   Psym(1)='kinss'
	   Psym(2)='kout'
	   Psym(3)='Sgmax'
	   Psym(4)='SC50g'
	   Psym(5)='Sbmax'
       Psym(6)='SC50b'
	   Psym(7)='kd'
	   Psym(8)='ga1'
	   Psym(9)='ga2'
	   Psym(10)='psi'
	   Psym(11)='Igmax'
	   Psym(12)='IC50g'
	   Psym(13)='Ibmax'
       Psym(14)='IC50b'
	   Psym(15)='kt'
	   Psym(16)='psi2'
	   Psym(17)='kt2'




CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
       PVsym(1)='int'
	   PVsym(2)='sig'


CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}     C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine DIFFEQ(T,X,XP)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 T,X(MaxNDE),XP(MaxNDE)
		Real*8 kout,Sgmax,SC50g,Sbmax,SC50b,kd,kinss,ga1,ga2,psi
		Real*8 gf1,gf2,gf3,bf1,bf2,bf3,Igmax,IC50g,Ibmax,IC50b
		Real*8 GF4,GF5,GF6,BF4,BF5,BF6,kt,PSI2,kt2
CC
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
		 kinss=P(1)
		 kout=P(2)
		 Sgmax=P(3)
		 SC50g=P(4)
		 Sbmax=P(5)
		 SC50b=P(6)
		 kd=P(7)
		 ga1=P(8)
		 ga2=P(9)
		 psi=P(10)
		 Igmax=P(11)
		 IC50g=P(12)
		 Ibmax=P(13)
		 IC50b=P(14)
		 kt=P(15)
		 psi2=P(16)
		 kt2=P(16)

		 
		 !NON-COMPETITIVE
		 gf1 = (psi*SC50g)
		 gf2 = psi2*IC50g !inhibit
 
		 bf1 = (psi*SC50b)
		 bf2 = psi2*IC50b!inhibit


		 XP(1)=X(5)-kout*X(1)
		 XP(2)=X(5)*(1-Igmax*6**ga1/(gf2**ga1+6**ga1))*
     c      (1-Ibmax*50**ga2/(bf2**ga2+50**ga2))-
     c      kout*(1+X(8))*
     c      (1+X(11))*X(2)
		 XP(3)=X(5)*(1-Igmax*108*ga1/(gf2**ga1+10**ga1))*
     c      (1-Ibmax*200**ga2/(bf2**ga2+200**ga2))-
     c      kout*(1+X(14))*
     c      (1+X(17))*X(3)
		 XP(4)=X(5)*(1-Igmax*20**ga1/(gf2**ga1+20**ga1))*
     c      (1-Ibmax*500**ga2/(bf2**ga2+500**ga2))-
     c      kout*(1+X(20))*
     c      (1+X(23))*X(4)
		
		 XP(5)=kd*X(5)*(1-X(5)/kinss)
		 !sgmax
		 XP(6)=kt*(Sgmax*6**ga1/(gf1**ga1+6**ga1)-X(6))
		 XP(7)=kt*(X(6)-X(7))
		 XP(8)=kt*(X(7)-X(8))
		 XP(9)=kt2*(Sbmax*50**ga2/(bf1**ga2+50**ga2)-X(9))
		 XP(10)=kt2*(X(9)-X(10))
		 XP(11)=kt*(X(10)-X(11))
		 XP(12)=kt*(Sgmax*10**ga1/(gf1**ga1+10**ga1)-X(12))
		 XP(13)=kt*(X(12)-X(13))
		 XP(14)=kt*(X(13)-X(14))
		 XP(15)=kt2*(Sbmax*200**ga2/(bf1**ga2+200**ga2)-X(15))
		 XP(16)=kt2*(X(15)-X(16))
		 XP(17)=kt2*(X(16)-X(17))
		 XP(18)=kt*(Sgmax*20**ga1/(gf1**ga1+20**ga1)-X(18))
		 XP(19)=kt*(X(18)-X(19))
		 XP(20)=kt*(X(19)-X(20))
		 XP(21)=kt2*(Sbmax*500**ga2/(bf1**ga2+500**ga2)-X(21))
		 XP(22)=kt2*(X(21)-X(22))
		 XP(23)=kt2*(X(22)-X(23))
		 
		 
		 
		 
		 
		 
		 
		 
		 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine OUTPUT(Y,T,X)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 Y(MaxNOE),T,X(MaxNDE),V


CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C
		 Y(1)=X(1)
		 Y(2)=X(2)
		 Y(3)=X(3)
		 Y(4)=X(4)

		 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine VARMOD(V,T,X,Y)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 V(MaxNOE),T,X(MaxNDE),Y(MaxNOE)

CC
C----------------------------------------------------------------------C
C   Enter Variance Model Equations Below                               C
C         {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                       C
C----c-----------------------------------------------------------------C

         V(1) = (PV(1) + PV(2)*Y(1))**2
		 V(2) = (PV(1) + PV(2)*Y(2))**2
		 V(3) = (PV(1) + PV(2)*Y(3))**2
		 V(4) = (PV(1) + PV(2)*Y(4))**2

		 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine COVMOD(Pmean, ICmean, PC)
C  Defines any covariate model equations (MLEM, ITS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 PC(MaxNCP)
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)

CC
C----------------------------------------------------------------------C
C     Enter # of Covariate Parameters                                  C
C----c-----------------------------------------------------------------C

        NCparam = 0    ! Enter # of Covariate Parameters.

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Covariate Params {eg: PCsym(1)='CLRenal'}         C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   For the Model Params. that Depend on Covariates Enter the Equation C
C         {e.g. Pmean(1) =  PC(1)*R(2) }                               C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine POPINIT(PmeanI,ICmeanI,PcovI,ICcovI, PCI)
C  Initial parameter values for population program parameters (ITS, MLEM)

        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 PmeanI(MaxNSP+MaxNDE), ICmeanI(MaxNDE)
        Real*8 PcovI(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcovI(MaxNDE,MaxNDE)
        Real*8 PCI(MaxNCP)

CC
C----------------------------------------------------------------------C
C  Enter Initial Values for Population Means                           C
C          {  e.g. PmeanI(1) = 10.0    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Initial Values for Pop. Covariance Matrix (Lower Triang.)    C
C         {  e.g. PcovI(2,1) = 0.25    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Values for Covariate Model Parameters                        C
C         {  e.g. PCI(1) = 2.0    }                                    C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
C  Parameter mean and covariance values for MAP estimation (ID,NPD,STS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
        Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(1) = 10.0    }                                C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Nonzero Elements of Covariance Matrix (Lower Triang.)       C
C         {  e.g. Pcov(2,1) = 0.25    }                                C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine SPARAM(PS,P,IC)
        Implicit None

        Include 'globals.inc'

        Real*8 PS(MaxNSECP), P(MaxNSP+MaxNDE), IC(MaxNDE) 

CC
C----------------------------------------------------------------------C
C   Enter Equations Defining Secondary Paramters                       C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C

        
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
	  	  
C######################################################################C

        Subroutine AMAT(A)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 A(MaxNDE,MaxNDE)

        DO I=1,Ndeqs
           Do J=1,Ndeqs
              A(I,J)=0.0D0
           End Do
        End Do

CC
C----------------------------------------------------------------------C
C   Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }   C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C