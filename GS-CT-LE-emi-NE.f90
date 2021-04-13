      program dsdrv1 

!     Pertains to ARPACK and P_ARPACK
!     Copyright (c) 1996-2008 Rice University.  
!     Developed by D.C. Sorensen, R.B. Lehoucq, C. Yang, and K. Maschhoff.
!     All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     The part related to optical spectra
!     is written by Xian-Kai Chen 
!     (chenxiankai2009@gmail.com)
!-----------------------------------------------------------------------
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer              maxn, maxnev, maxncv, ldv
      parameter            (maxn=2*20*200, maxnev=1000, maxncv=2000, &
     &                     ldv=maxn )
      integer,parameter :: mm=20
      integer,parameter :: kk=200
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      Double precision  &
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),&
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),&
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, j,& 
     &                 nx, nconv, maxitr, mode, ishfts
      logical          rvec
      Double precision   &   
     &                 tol, sigma
      integer           jj, je1, je2
      integer           k1, k2, k3, k4 
      character(3)      c  
      real(8)           e2, ee2,  av_abs, av_emi 
      integer,parameter :: maxe2=2001
      
      real(8)           A_e2(maxe2),  t_av_emi(maxe2)
      real(8)           r2     
      
      real(8)           Ttdm
      real(8)           tdmle, tdmct                                                  ! tdm: transition dipole moment for diabatic states                                   
     
      real(8)           BolzES, BolzGS
      real(8)           freqm, freqk
      real(8)           freqE(mm*kk)
      real(8)           kbT                                                           ! Bolzmann constant * temparature                                          

!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision  &
     &                 zero
      parameter        (zero = 0.0D+0)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  &         
     &                 dnrm2
      external         dnrm2, daxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic        abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number NX is the number of interior points     |
!     | in the discretization of the 2-dimensional         |
!     | Laplacian on the unit square with zero Dirichlet   |
!     | boundary condition.  The number N(=NX*NX) is the   |
!     | dimension of the matrix.  A standard eigenvalue    |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
!
      n = 2*20*200
      nev =  1000
      ncv =  2000
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'SA'
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 1
!      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

      jj=4                                                                             
      write(c,'(i3)') jj
      open(12,file='3_states_energy_'//trim(adjustl(c))//'.dat',status='unknown')

 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,  &
     &                 lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            call av (n,  workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if 
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, &
     &        iparam, ipntr, workd, workl, lworkl, ierr )
     
      write(12,15) d(1:nev,1)
  15  format(f18.9) 
         
      open(14,file='pbt7-th-emi_'//trim(adjustl(c))//'.dat',status='unknown')                      !spectral function
      
      kbT   = 200.0d0
      freqm = 1200.0d0
      freqk = 100.0d0

      tdmle = 20.0d0
      tdmct = 0.0d0

      r2  = 100.0d0
      ee2 = 10.0d0

      BolzES = 0.0d0
      do k1 = 1, nev
       BolzES = BolzES + exp(-d(k1,1)/kbT)
      enddo

      freqE  = 0.0d0
      BolzGS = 0.0d0
      do k2 = 0, mm-1
        do k3= 0, kk-1    
          freqE(k2*kk+k3+1) =  freqm*(0.5+k2) + freqk*(0.5+k3)        
          BolzGS = BolzGS + exp(-(freqm*(0.5+k2) + freqk*(0.5+k3))/kbT) 
        enddo
      enddo

      A_e2      =  0.0d0
      t_av_emi  =  0.0d0

         e2 = -ee2
         do je2=1,maxe2

           e2 = e2 + ee2
           av_emi  = 0.0d0
!=========================================emi and abs======================================================================
          do k1 = 1,nev                                                     !ES
           do k2 = 1,mm                                                    !GS, vibration states
            do k3 = 1,kk
              if((d(k1,1)-(freqm*(0.5+k2-1) + freqk*(0.5+k3-1))).gt.0.0d0)    then
                    
                     Ttdm  = 0.0d0
                     Ttdm  = 0.5*(tdmle*v((k2-1)*kk+k3,k1))**2 + 0.5*(tdmct*v(mm*kk+(k2-1)*kk+k3,k1))**2
           
                    av_emi  = av_emi  + (e2**3) * (exp(-d(k1,1)/kbT)/BolzES) * (Ttdm) * exp(-(e2-(d(k1,1)-(freqm*(0.5+k2-1) + freqk*(0.5+k3-1))))**2/r2**2)    
              else
                     goto 617 
              endif
            enddo
617         continue
           enddo 
          enddo                                
!=========================================================================================================================
         A_e2(je2)      = e2
         t_av_emi(je2)  = av_emi           
         enddo

         t_av_emi  = t_av_emi / maxval(t_av_emi)

         write(14,17) (A_e2(je2), t_av_emi(je2),  je2=1,maxe2)
17       format(2f18.9)

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                call av(n, v(1,j), ax)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
!
 20          continue
!
!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
!
             call dmout(6, nconv, 2, d, maxncv, -6,  &
     &            'Ritz values and relative residuals')
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit  &
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      
!
         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',   &
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue
!
      end
!-------------------------------------------------------------------
SUBROUTINE av(n,v,w)
    implicit none
    integer n
    integer,parameter   :: maxm=19, maxk=199
    double precision v(n), w(n)
    call multv(maxk,maxm,v,w)
end SUBROUTINE av

SUBROUTINE multv(maxk,maxm,v,w)
    IMPLICIT NONE
    INTEGER, parameter              ::  maxn=2
    INTEGER, intent(in)             ::  maxm, maxk
    DOUBLE PRECISION, intent(in)    ::  v(0:maxk,0:maxm,maxn)
    DOUBLE PRECISION, intent(out)   ::  w(0:maxk,0:maxm,maxn)
    ! Defining constants of vibrational Hamiltonian
    REAL(8)                             tctle                                  ! transfer integral             
    REAL(8)                             g1, g2, g3, g4                         ! e-v coupling           
    REAL(8)                             f1, f2                                 ! freq              
    REAL(8)                             AEle, AEct                             ! AE2: diabtic S1 energy; AE3: diabatic CT energy                                                          


    ! quantum numbers for rows of H
    INTEGER m, k

    AEle  = 13389.0d0       ! / cm-1                        
    AEct  = 13389.0d0       ! / cm-1    3333 ~ 6333, 500
   
    tctle = 0.0d0        
            
    g1  = -0.75d0
    g2  = -2.3d0
    g3  = 0.0d0
    g4  = 0.0d0

    f1  = 1200.0d0        ! / cm-1            
    f2  = 100.0d0        ! / cm-1          

    w = 0d0
      
    do m=0,maxm
     do k=0,maxk
        
        
                    w(k,m,1) = w(k,m,1) + v(k,m,1)*((m+0.5)*f1+(k+0.5)*f2+AEle+g1*g1*f1+g2*g2*f2)
        if(m>0)     w(k,m,1) = w(k,m,1) + v(k,m-1,1)*g1*f1*sqrt(dble(m))
        if(m<maxm)  w(k,m,1) = w(k,m,1) + v(k,m+1,1)*g1*f1*sqrt(dble(m+1))
        if(k>0)     w(k,m,1) = w(k,m,1) + v(k-1,m,1)*g2*f2*sqrt(dble(k))
        if(k<maxk)  w(k,m,1) = w(k,m,1) + v(k+1,m,1)*g2*f2*sqrt(dble(k+1))
        
                    w(k,m,2) = w(k,m,2) + v(k,m,2)*((m+0.5)*f1+(k+0.5)*f2+AEct+g3*g3*f1+g4*g4*f2)
        if(m>0)     w(k,m,2) = w(k,m,2) + v(k,m-1,2)*g3*f1*sqrt(dble(m))
        if(m<maxm)  w(k,m,2) = w(k,m,2) + v(k,m+1,2)*g3*f1*sqrt(dble(m+1))                                                 
        if(k>0)     w(k,m,2) = w(k,m,2) + v(k-1,m,2)*g4*f2*sqrt(dble(k))
        if(k<maxk)  w(k,m,2) = w(k,m,2) + v(k+1,m,2)*g4*f2*sqrt(dble(k+1))                                   
           
                    w(k,m,1) = w(k,m,1) + v(k,m,2)*tctle   
                    w(k,m,2) = w(k,m,2) + v(k,m,1)*tctle            
                    
     enddo 
    enddo 
END SUBROUTINE multv


