!
!>\file mp_ufs.F90
!! This file is the ufs Micorphysics shceme.
!! author Songyou Hong, songyou.hong@noaa.gov
  module mp_ufs
   use module_mp_radar
   use machine , only : kind_phys
   implicit none
   public :: mp_ufs_init, mp_ufs_run, mp_ufs_finalize
   real(kind=kind_phys), parameter, private :: &
          dtcldcr = 180. ,& !< maximum time step for sub-cycling loops [s]
          n0r = 8.e6 ,& !< intercept parameter for rain [m-4]
          n0s = 2.e6 ,& !< intercept parameter constant for snow [m-4]
          n0g = 4.e6 ,& !< intercept parameter graupel
          n0smax = 1.e11 ,& !< maximum n0s (t=-90c unlimited)
          alpha = .12 ,& !< .122 exponen factor for n0s computation
          avtr = 841.9 ,& !< a constant for terminal velocity of rain
          bvtr = 0.8 ,& !< a constant for terminal velocity of rain
          avts = 11.72 ,& !< a constant for terminal velocity of snow
          bvts = .41 ,& !< a constant for terminal velocity of snow
          avtg = 330. ,& !< a constant for terminal velocity of graupel
          bvtg = 0.8 ,& !< a constant for terminal velocity of graupel
          deng = 500. ,& !< density of graupel
          lamdasmax = 1.e5 ,& !< maximum value of slope parameter for snow
          lamdagmax = 6.e4 ,& !< maximum value of slope parameter for graupel
          r0 = .8e-5 ,& !< critical droplet radius in autoconversion [m]
          peaut = .55 ,& !< collection efficiency in autoconversion
          xncr = 3.e8 ,& !< number concentration of cloud droplets [m-3]
          xmyu = 1.718e-5,& !< the dynamic viscosity [kgm-1s-1]
          dicon = 11.9 ,& !< constant for the cloud-ice diamter
          dimax = 500.e-6 ,& !< maximum value for the cloud-ice diamter [m]
          pfrz1 = 100. ,& !< constant in biggs freezing
          pfrz2 = 0.66 ,& !< constant in biggs freezing
          qcmin = 1.e-12 ,& !< minimun values for qr, qs, and qg [kgkg-1]
          qrmin = 1.e-9 ,& !< minimun values for qr, qs, and qg [kgkg-1]
          eacrc = 1.0 ,& !< snow/cloud-water collection efficiency
          dens = 100.0 ,& !< density of snow
          qs0 = 6.e-4 !< threshold amount for aggretion to occur
   real(kind=kind_phys), parameter, private :: &
          drcoeff = (24.)**(.3333333) ,& !< factor for raindrop diameter
          lamdarmax = 5.e4 ,& !< maximum value of slope parameter for rain
          lamdarmin = 2.e3 ,& !< minimum value of slope parameter for rain
          lamdacmax = 5.e5 ,& !< maximum value of slope parameter for cloud
          lamdacmin = 2.e4 ,& !< minimum value of slope parameter for cloud
          ncmin = 1.e1 ,& !<
          nrmin = 1.e-2 ,& !<
          satmax = 0.0048 ,& !<
          actk = 0.6 ,& !<
          actr = 1.5 ,& !<
          ncrk1 = 3.03e3 ,& !<
          ncrk2 = 2.59e15 ,& !<
          ncr2micro =2.39e14 ,& !<
          ncr82micro=3.47e9 ,& !<
          ncr5mm = 1.53e4 ,& !<
          di5000 = 5000.e-6 ,& !<
          di2000 = 2000.e-6 ,& !<
          di600 = 6.e-4 ,& !<
          di100 = 1.e-4 ,& !<
          di82 = 82.e-6 ,& !<
          di15 = 15.e-6 ,& !<
          di2 = 2.e-6 !<
   real(kind=kind_phys), parameter, private :: &
          dust0 = 2.5 ,& !< background dust (microgram m-3)
          ccn_0 = 50.e+6 ,& !< background ccn at t = 0 (number m-3)
          ccnmax = 20000.e+6,& !< maximum ccn (number m-3)
          ccnmin = 20.e+6 !< minimum ccn (number m-3)
   real(kind=kind_phys), parameter, private :: &
          qmin = 1.e-15 ,& !< minimun values for tracers [kgkg-1]
          xls = 2.8340e+6 ,& !< latent hear for sublimation
          xlv0 = 2.5010e+6 ,&
          xlf0 = 3.3370e+5 ,&
          den0 = 1.28 ,&
          denr = 1000. ,&
          cliq = 4.1900e+3 ,&
          cice = 2.1060e+3 ,&
          psat = 6.1078e+2
   real(kind=kind_phys), parameter, private :: &
          recmin = 2.51e-6 ,&
          recmax = 50.e-6 ,&
          reimin = 5.01e-6 ,&
          reimax = 125.e-6 ,&
          resmin = 25.e-6 ,&
          resmax = 999.e-6
!
   real(kind=kind_phys), parameter, private :: &
          sedi_semi_cfl = 10.0 ,&
          zero_0 = 0.0 ,&
          one_1 = 1.0
   logical, parameter, private :: flgzero = .true.
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   real(kind=kind_phys), save :: &
             qc_ocean,qc_land,qck1,pidnc,bvtr1,bvtr2,bvtr3,bvtr4,bvtr5, &
             bvtr6,bvtr7, bvtr2o5,bvtr3o5, &
             g1pbr,g2pbr,g3pbr,g4pbr,g5pbr,g6pbr,g7pbr, &
             g5pbro2,g7pbro2,pi,pisq, &
             pvtr,pvtrn,eacrr,pacrr,pidn0r,pidnr, &
             precr1,precr2,xmmax,roqimax,bvts1,bvts2, &
             bvts3,bvts4,g1pbs,g3pbs,g4pbs,g5pbso2, &
             pvts,pacrs,precs1,precs2,pidn0s,xlv1,pacrc, &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,g3pbg,g4pbg, &
             g5pbgo2,pvtg,pacrg,precg1,precg2,pidn0g, &
             rslopecmax,rslopec2max,rslopec3max, &
             rslopermax,rslopesmax,rslopegmax, &
             rsloperbmax,rslopesbmax,rslopegbmax, &
             rsloper2max,rslopes2max,rslopeg2max, &
             rsloper3max,rslopes3max,rslopeg3max
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!! \section arg_table_mp_ufs_init Argument Table
!! \htmlinclude mp_ufs_init.html
!!
!-------------------------------------------------------------------------------
      subroutine mp_ufs_init (ncol, nlev, qc, qr, qi, qs, qg, &
                 nn, nc, nr, mpirank, mpiroot, errmsg, errflg)
         implicit none
         ! Interface variables
         integer, intent(in ) :: ncol
         integer, intent(in ) :: nlev
         real(kind_phys), intent(inout) :: qc(ncol,nlev)
         real(kind_phys), intent(inout) :: qr(ncol,nlev)
         real(kind_phys), intent(inout) :: qi(ncol,nlev)
         real(kind_phys), intent(inout) :: qs(ncol,nlev)
         real(kind_phys), intent(inout) :: qg(ncol,nlev)
         real(kind_phys), intent(inout) :: nn(ncol,nlev)
         real(kind_phys), intent(inout) :: nc(ncol,nlev)
         real(kind_phys), intent(inout) :: nr(ncol,nlev)
         integer, intent(in) :: mpirank
         integer, intent(in) :: mpiroot
         character(len=*), intent( out) :: errmsg
         integer, intent( out) :: errflg
         real(kind_phys), parameter :: ccn0 = 100.e+6
         real(kind_phys), parameter :: den0 = 1.28
         real(kind_phys), parameter :: denr = 1000.
         real(kind_phys), parameter :: dens = 100.
         real(kind_phys), parameter :: cl = 4.1855e+3
         real(kind_phys), parameter :: cpv = 1.8460e+3
!
         integer :: i, k
         ! Initialize the ccpp error handling variables
         errmsg = ''
         errflg = 0
         ! *DH temporary
         if (mpirank==mpiroot) then
            write(0,*) ' -------------------------------------------------------------------'
            write(0,*) ' --- WARNING --- the CCPP UFS Microphysics scheme is currently under development, use at your own risk --- WARNING ---'
            write(0,*) ' -------------------------------------------------------------------'
         end if
         ! *DH temporary
!
         qc(:,:) = 0.0 ; qr(:,:) = 0.0
         qi(:,:) = 0.0 ; qs(:,:) = 0.0 ; qg(:,:) = 0.0
!
         nn(:,:) = ccn0
         nc(:,:) = 0.0 ; nr(:,:) = 0.0
!
         call ufsmpinit(den0,denr,dens,cl,cpv)
!
      end subroutine mp_ufs_init
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
      subroutine mp_ufs_finalize()
!-------------------------------------------------------------------------------
      end subroutine mp_ufs_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!> \section arg_table_mp_ufs_run Argument Table
!! \htmlinclude mp_ufs_run.html
!!
subroutine mp_ufs_run(ncol, nlev, t1, q1 &
                         ,qc1, qr1, qi1, qs1, qg1 &
                         ,nn1,nc1,nr1 &
                         ,p1, phii, delt, g, cpd, cpv, rd, rv, t0c &
                         ,ep1, ep2 &
                         ,area &
                         ,slmsk &
                         ,hpbl &
                         ,rain,snow, graupel &
! ,qmin, xls, xlv0, xlf0, den0, denr &
! ,cliq,cice,psat &
! ,lat &
! ,sr, snow, snowncv &
! ,graupel,graupelncv &
                         ,refl_10cm, re_cloud, re_ice, re_snow &
                         ,fhswr, kdt, aerfld &
                         ,errmsg, errflg &
                                                                            )
!
!===================================================================
!! Unified forecast system microphyiscs scheme (ufs mp) utilizes the
!! in-cloud microphysics concept for production terms (kim and hong 2018).
!! usm mp assumes that ice nuclei number concentration is a function of
!! temperature, and seperate assumption is developed, in which ice crystal
!! number concentration is a function of ice amount (hong et al. 2004).
!! double moment microphysics for warm microphysics is adapted from lim and
!! hong (2010), having prognostic number concentrations for cloud, rain, and
!! cloud condensation nuclei (ccn) concentrations. semi-lagrangian
!! sedimentation of juang and hong (2010) is also re-configured for
!! computational efficiency and numerical accuracy. all productions terms
!! are optimized by introducing cloud top definition for hydrometeors.
!!
!! all units are in m.k.s. and source/sink terms in kgkg-1s-1.
!!
!! ufs mp scheme, coded by songyou hong (noaa/esrl/psl)
!! summer, fall and winter 2021
!! wdm mp scheme, coded by kyo-sun lim and songyou hong (yonsei univ.)
!! fall 2008
!! wsm mp scheme, coded by songyou hong, jeong-ock jade lim (yonsei univ.),
!! jimy dudhia, and shu-hua chen (ncar)
!! summer 2002
!! ncep mp scheme, coded by songyou hong, henry juang, qingyun zhao (ncep),
!! summer 1997
!! references :
!! kim and and hong (kh, 2018) j. atmos. sci.
!! juang and hong (jh, 2010) mon. wea. rev.
!! lim and hong (lh, 2010) mon. wea. rev.
!! dudhia, hong and lim (dhl, 2008) j. meteor. soc. japan
!! hong and lim (hl, 2006) j. korean meteor. soc.
!! hong, dudhia, and chen (hdc, 2004) mon. wea. rev.
!! hong, juang, and zhao (hjz, 1998) mon. wea. rev.
!!
!! \section structure
!! --- ufs_mps --- |-slope_rain slope_snow slope_graupel
!! |
!! |- semi_lagrange_sedim
!! |
!! |- cldf_mps_diag
!! |
!! |- ufs_mp_reflectivity
!! |
!! |- effective_radius
!!
!! \ section input and output
!! input :
!! delt - timestep
!! g, cpd, cpv, t0c, den0, - constant
!! rd, rv, ep1, ep2, qmin,
!! xls, xlv0, xlf0, denr,
!! cliq, cice, psat
!! ids, ide, jds, jde, kds, kde - dimension
!! ims, ime, jms, jme, kms, kme
!! its, ite, jts, jte, kts, kte
!! ncloud - number of hydrometeor
!! p - pressure, pa
!! delz - depth of model layer, m
!!
!! inout :
!! t1 - temperautre
!! q1 - specific humidity
!! q2 - mixing ratio of cloud, rain, ice, and snow
!! qc, qr, qi, qc
!! rain, rainncv - precipitation
!! sr - ratio of snow to rain
!!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer , intent(in ) :: ncol, nlev, kdt
   real(kind=kind_phys) , intent(in ) :: fhswr
   real(kind=kind_phys) , intent(in ) :: delt
   real(kind=kind_phys) , intent(in ) :: g, cpd, cpv, t0c, &
                                                             rd, rv, ep1, ep2
   integer, dimension( ncol ) , intent(in ) :: slmsk
   real(kind=kind_phys), dimension( ncol ), intent(in ) :: hpbl
   real(kind=kind_phys), dimension( ncol,nlev ), intent(in ) :: p1
   real(kind=kind_phys), dimension( ncol,nlev + 1), intent(in ) :: phii
   real(kind=kind_phys), dimension( ncol,nlev ), intent(inout) :: q1
   real(kind=kind_phys), dimension( ncol,nlev ), intent(inout) :: t1
   real(kind=kind_phys), dimension( ncol ) , intent(in ) :: area
   real(kind=kind_phys), dimension( ncol,nlev ), intent(inout) :: qc1, qr1, qi1, qs1, qg1
   real(kind=kind_phys), dimension( ncol,nlev ) , intent(inout) :: nn1, nc1, nr1
   real(kind=kind_phys), dimension( ncol ) , intent(inout) :: rain
   real(kind=kind_phys), dimension( ncol ), optional , intent(inout) :: snow
   real(kind=kind_phys), dimension( ncol ), optional , intent(inout) :: graupel
   real(kind=kind_phys), dimension( ncol , nlev ), intent(inout) :: &
                               refl_10cm, re_cloud, re_ice, re_snow
   real(kind=kind_phys), intent(in ) :: aerfld(:,:,:)
   character(len=*), intent(out) :: errmsg
   integer, intent(out) :: errflg
!
! local variables
!
   real(kind=kind_phys), dimension( ncol , nlev ) :: delz1, tv
   real(kind=kind_phys), dimension( ncol ) :: dxmeter
   real(kind=kind_phys), dimension( ncol ) :: cbot
   real(kind=kind_phys), dimension( nlev , ncol ) :: &
                                                                q, &
                                                              den, &
                                                                t, &
                                                                p, &
                                                             delz, &
                                                             dend, &
                                                               zl, &
                                                           denfac, &
                                                              xlv, &
                                                              xlf, &
                                                              cpm
   real(kind=kind_phys), dimension( nlev , ncol ) :: &
                                                              qc0
   real(kind=kind_phys), dimension( nlev , ncol , 2) :: &
                                                              qci
   real(kind=kind_phys), dimension( nlev , ncol , 3) :: &
                                                              qrs
   real(kind=kind_phys), dimension( nlev, 3) :: &
                                                           rslope, &
                                                          rslope2, &
                                                          rslope3, &
                                                          rslopeb
   real(kind=kind_phys), dimension( nlev , ncol , 3) :: ncr
   real(kind=kind_phys), dimension( nlev ) :: &
                                                           rh_mul, &
                                                           rh_ice, &
                                                         qsat_mul, &
                                                         qsat_ice, &
                                                              vtr, &
                                                              vts, &
                                                              vtg, &
                                                             vtsg, &
                                                              vti, &
                                                               ni, &
                                                            sumsg, &
                                                            denqr, &
                                                            denqs, &
                                                            denqg, &
                                                            denqi, &
                                                           n0sfac, &
                                                           ab_mul, &
                                                           ab_ice, &
                                                             diai, &
                                                           venfac
   real(kind=kind_phys), dimension( nlev ) :: &
                                                          rslopec, &
                                                         rslopec2, &
                                                         rslopec3, &
                                                             diac, &
                                                             diar, &
                                                            denqn
   real(kind=kind_phys), dimension( nlev ) :: &
                                                            pigen, &
                                                            pcond, &
                                                            prevp, &
                                                            psevp, &
                                                            pgevp, &
                                                            pidep, &
                                                            psdep, &
                                                            pgdep, &
                                                            praut, &
                                                            psaut, &
                                                            pgaut, &
                                                            piacr, &
                                                            pracw, &
                                                            praci, &
                                                            pracs, &
                                                            psacw, &
                                                            psaci, &
                                                            psacr, &
                                                            pgacw, &
                                                            pgaci, &
                                                            pgacr, &
                                                            pgacs, &
                                                            paacw, &
                                                            psmlt, &
                                                            pgmlt, &
                                                            pseml, &
                                                            pgeml
   real(kind=kind_phys), dimension( nlev ) :: &
                                                            pcact, &
                                                            nraut, &
                                                            nracw, &
                                                            ncevp, &
                                                            nccol, &
                                                            nrcol, &
                                                            nsacw, &
                                                            ngacw, &
                                                            niacr, &
                                                            nsacr, &
                                                            ngacr, &
                                                            naacw, &
                                                            nseml, &
                                                            ngeml, &
                                                            ncact
   real(kind=kind_phys), dimension( nlev ) :: &
                                                           satrdt, &
                                                           supice, &
                                                           supsat, &
                                                           tcelci, &
                                                             cldf
   real(kind=kind_phys) :: &
                                                           qrpath, &
                                                           qspath, &
                                                           qgpath, &
                                                           qipath, &
                                                         precip_r, &
                                                         precip_s, &
                                                         precip_g, &
                                                         precip_i
   logical, dimension( nlev ) :: &
                                                            ifsat, &
                                                            ifice
!
! for reflectivity....
!
   real(kind=kind_phys), dimension( nlev ) :: &
                                                              dbz, &
                                                            re_qc, &
                                                            re_qi, &
                                                            re_qs
!
! miscellaneous variables
!
   real(kind=kind_phys) :: &
             rdtcld,conden, x, y, z, a, b, c, d, e, &
             coeres, dtcld, massi, eacrs, acrfac, egi, &
             qimax, ni0, roqi0, &
             precip_sum, precip_ice, factor, source, total, &
             value, pfrzdtc, pfrzdtr, &
             tstepsnow, tstepgraupel, xal, xbl, alpha2, &
             delta2, delta3, dldti, xb, xai, tr, xbi, xa, hvap, &
             cvap, hsub, dldt, ttp, logtr, dtcfl, temp
   real(kind=kind_phys) :: gfac, sfac, nfrzdtr, nfrzdtc, qnpath
   real(kind=kind_phys), dimension( nlev ) :: qrconcr, qrcon, taucon
   real(kind=kind_phys), dimension( nlev ) :: vtn
   integer :: nstep, niter, lond, latd, &
              i, j, k, loop, loops, n, idim, kdim
   integer :: ktopini, ktopmax, ktop
   integer :: ktopqc, ktopqi, ktopqr, ktopqs, ktopqg, ktoprh
   integer :: its, ite, kts, kte
   integer :: lat
   logical :: flgcld, lqc, lqi, lqr, lqs, lqg
!
!-------------------------------------------------------------------------------
! preparation ....
!
   its = 1
   ite = ncol
   kts = 1
   kte = nlev
!--- initialize CCPP error handling variables
   errmsg = ''
   errflg = 0
   idim = ite-its+1
   kdim = kte-kts+1
   ktopini = kte - 1
   lond = (ite-its)/2 + 1
   lat = 1
   latd = 1
   dxmeter(:) = sqrt(area(:))
!-------------------------------------------------------------------------------
! assign passing variables to local arrays
!
   zl(:,:) = 0.; cbot(:) = 0.
   do k = kts, kte
     do i = its, ite
       tv(i,k) = t1(i,k)*(1.0+ep1*q1(i,k))
     enddo
   enddo
!
   delz1(:,:) = 0.0
   do k = kts, kte
     do i = its, ite
       delz1(i,k) = (phii(i,k+1)-phii(i,k))/g
     enddo
   enddo
   do k = kts, kte
     do i = its, ite
       t(k,i) = t1(i,k)
       q(k,i) = q1(i,k)
       p(k,i) = p1(i,k)
       delz(k,i) = delz1(i,k)
       den(k,i) = p(k,i)/(rd*tv(i,k))
       dend(k,i) = (p(k,i)/t(k,i)-den(k,i)*rv)/(rd-rv) ! dry density
       denfac(k,i) = sqrt(den0/den(k,i))
     enddo
   enddo
   zl(1,:) = delz(1,:)*0.5
   do k = kts+1, kte
     do i = its, ite
       zl(k,i) = zl(k-1,i)+(delz(k-1,i)+delz(k,i))*0.5
     enddo
   enddo
   do i = its, ite
     cbot(i) = max(hpbl(i),500.)
   enddo
!
!-------- initlize CCN -----
!
   if (kdt==1) then
     do k = kts, kte
       do i = its, ite
         if(slmsk(i)==0) then
! oceans nn1(i,k)=10.**(2.41+0.5*log(aerfld(SO4))+0.13*log(1.3*aerfld(OCphilic))&
! +0.05*log(aerfld(Seasalt)))
           nn1(i,k) = 10.**(2.41+0.5*log(aerfld(i,k,11)*den(k,i)*1.e9) &
                              +0.13*log(1.3*aerfld(i,k,15)*den(k,i)*1.e9) &
                              +0.05*log(aerfld(i,k,6)*den(k,i)*1.e9 &
                                       +aerfld(i,k,7)*den(k,i)*1.e9 &
                                       +aerfld(i,k,8)*den(k,i)*1.e9))*1.e6
         else
! land nn1(i,k)=10.**(2.41+0.5*log(aerfld(SO4))+0.13*log(1.3*aerfld(OCphilic)))
           nn1(i,k) = 10.**(2.41+0.5*log(aerfld(i,k,11)*den(k,i)*1.e9) &
                    +0.13*log(1.3*aerfld(i,k,15)*den(k,i)*1.e9))*1.e6
         endif
         nn1(i,k) = nn1(i,k) + ccn_0 ! Add a background value for CCN
       enddo
     enddo
   endif
!-------------------------------------------------------------------------------
! padding 0 for negative values generated by dynamics - optional
!
   ktop = ktopini
   do k = kts,ktop
     do i = its,ite
       if(flgzero) then
         qci(k,i,1) = max(qc1(i,k),0.0)
         qrs(k,i,1) = max(qr1(i,k),0.0)
         qci(k,i,2) = max(qi1(i,k),0.0)
         qrs(k,i,2) = max(qs1(i,k),0.0)
         qrs(k,i,3) = max(qg1(i,k),0.0)
         ncr(k,i,1) = min(max(nn1(i,k),ccnmin),ccnmax)
         ncr(k,i,2) = max(nc1(i,k),0.0)
         ncr(k,i,3) = max(nr1(i,k),0.0)
       else
         qci(k,i,1) = qc1(i,k)
         qrs(k,i,1) = qr1(i,k)
         qci(k,i,2) = qi1(i,k)
         qrs(k,i,2) = qs1(i,k)
         qrs(k,i,3) = qg1(i,k)
         ncr(k,i,1) = min(max(nn1(i,k),ccnmin),ccnmax)
         ncr(k,i,2) = nc1(i,k)
         ncr(k,i,3) = nr1(i,k)
       endif
     enddo
   enddo
!
   do i = its,ite
     do k = ktop+1,kte
       qci(k,i,1) = 0.0
       qrs(k,i,1) = 0.0
       qci(k,i,2) = 0.0
       qrs(k,i,2) = 0.0
       qrs(k,i,3) = 0.0
     enddo
   enddo
!
! initialize the surface rain, snow, graupel and etc....
!
!
   rain(:) = 0.; snow(:) = 0.; graupel(:) = 0.
   cpm(:,:) = 0.; xlv(:,:) = 0.; xlf(:,:) = 0.
!-------------------------------------------------------------------------------
! latent heat for phase changes and heat capacity. emanuel(1994)
!
   do k = kts, ktop
     do i = its, ite
       cpm(k,i) = cpd*(1.-max(q(k,i),qmin)) + max(q(k,i),qmin)*cpv
       xlv(k,i) = xlv0 - xlv1*(t(k,i)-t0c)
       xlf(k,i) = xls - xlv(k,i)
     enddo
   enddo
!
   hsub = xls
   hvap = xlv0
   cvap = cpv
   ttp = t0c+0.01
   dldt = cvap-cliq
   xa = -dldt/rv
   xb = xa+hvap/(rv*ttp)
   dldti= cvap-cice
   xai = -dldti/rv
   xbi = xai+hsub/(rv*ttp)
!===============================================================================
!
! inner loop wih the time step of dtcldcr (default = 180 sec)
!
!===============================================================================
! compute the sub-cycling time steps.
!
   loops = max(ceiling(delt/dtcldcr),1)
   dtcld = delt/loops
   if(delt<=dtcldcr) dtcld = delt
   rdtcld = 1./dtcld
!
 inner_loop : do loop = 1, loops
!
!-------------------------------------------------------------------------------
 i_loop : do i = its, ite ! i-loop for one-dimensional code
!
   qsat_mul(:) = 0.; rh_mul(:) = 0.
   qsat_ice(:) = 0.; rh_ice(:) = 0.
   ktopqc = 0; ktopqi = 0; ktopqr = 0; ktopqs = 0; ktopqg = 0; ktoprh = 0
!
   ktop = ktopini
   do k = kts, ktop
     if(t(k,i)<ttp) then
       xal = xai
       xbl = xbi
     else
       xal = xa
       xbl = xb
     endif
     tr=ttp/t(k,i)
     logtr=log(tr)
     qsat_mul(k) = psat*exp(logtr*(xa)+xb*(1.-tr))
     qsat_mul(k) = min(qsat_mul(k),0.99*p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k) = max(q(k,i) / qsat_mul(k),qmin)
     qsat_ice(k) = psat*exp(logtr*(xal)+xbl*(1.-tr))
     qsat_ice(k) = min(qsat_ice(k),0.99*p(k,i))
     qsat_ice(k) = ep2 * qsat_ice(k) / (p(k,i) - qsat_ice(k))
     qsat_ice(k) = max(qsat_ice(k),qmin)
     rh_ice(k) = max(q(k,i) / qsat_ice(k),qmin)
   enddo
!
!-------------------------------------------------------------------------------
! ab: the thermodynamic term in the denominator associated with
! heat conduction and vapor diffusion
! venfac: parameter associated with the ventilation effects
!
   do k = ktop, kts, -1
! ab_mul(k) = ab_mul(xlv(k,i),p(k,i),t(k,i),den(k,i),qsat_mul(k))
       ab_mul(k) = ((((den(k,i))*(xlv(k,i))*(xlv(k,i)))*((t(k,i))+120.) &
                    *(den(k,i))) /(1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i)))) &
                    *(den(k,i))* (rv*(t(k,i))*(t(k,i))))) &
                    + p(k,i)/((qsat_mul(k))*(8.794e-5*exp(log(t(k,i))*(1.81))))
! ab_ice(k) = ab_ice(xls,p(k,i),t(k,i),den(k,i),qsat_ice(k))
       ab_ice(k) =((((den(k,i))*(xls)*(xls))*((t(k,i))+120.) &
                    *(den(k,i))) /(1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i)))) &
                    *(den(k,i))* (rv*(t(k,i))*(t(k,i)))) &
                    + p(k,i)/(qsat_ice(k)*(8.794e-5*exp(log(t(k,i))*(1.81)))))
! venfac(k) = venfac(p(k,i),t(k,i),den(k,i))
       venfac(k) = (exp(.3333333*log(((1.496e-6 * ((t(k,i))*sqrt(t(k,i)))) &
                  *p(k,i))/(((t(k,i))+120.)*den(k,i)*(8.794e-5 &
                  *exp(log(t(k,i))*(1.81)))))) *sqrt(sqrt(den0/(den(k,i))))) &
                   /sqrt((1.496e-6*((t(k,i))*sqrt(t(k,i)))) &
                   /(((t(k,i))+120.)*den(k,i)))
   enddo
!-------------------------------------------------------------------------------
! re-define cloud top for numerical efficiencies
!
   lqc = .false.
   lqi = .false.
   lqr = .false.
   lqs = .false.
   lqg = .false.
   flgcld = .false.
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,qci(:,i,2),zero_0,ktopqi)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,2),zero_0,ktopqs)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,3),zero_0,ktopqg)
   call find_cloud_top(1,kdim,ktopini,rh_ice(:), one_1,ktoprh)
   if(ktopqc>0.0) lqc = .true.
   if(ktopqi>0.0) lqi = .true.
   if(ktopqr>0.0) lqr = .true.
   if(ktopqs>0.0) lqs = .true.
   if(ktopqg>0.0) lqg = .true.
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktoprh)
!-------------------------------------------------------------------------------
! early checkout
!
   if(lqc .or. lqi .or. lqr .or. lqs .or. lqg) flgcld = .true.
   if((.not.flgcld) .and. ktoprh>0.0) flgcld = .true.
   if(.not.flgcld) then
     cycle i_loop
   endif
!-------------------------------------------------------------------------------
! initialize the local variables
!
   temp = 0. ; value = 0. ; total = 0. ; source = 0.
   pcond(:) = 0. ; pigen(:) = 0.
   praut(:) = 0. ; psaut(:) = 0. ; pgaut(:) = 0.
   pidep(:) = 0. ; psdep(:) = 0. ; pgdep(:) = 0.
   pracw(:) = 0. ; praci(:) = 0. ; pracs(:) = 0.
   psacw(:) = 0. ; psaci(:) = 0. ; psacr(:) = 0.
   pgacw(:) = 0. ; pgaci(:) = 0. ; pgacr(:) = 0. ; pgacs(:) = 0.
   paacw(:) = 0. ; piacr(:) = 0.
   psmlt(:) = 0. ; pgmlt(:) = 0. ; pseml(:) = 0. ; pgeml(:) = 0.
   prevp(:) = 0. ; psevp(:) = 0. ; pgevp(:) = 0.
   denqr(:) = 0. ; denqs(:) = 0. ; denqi(:) = 0.
   vtr(:) = 0. ; vts(:) = 0. ; vtg(:) = 0. ; vti(:) = 0.
   qrpath = 0. ; qspath = 0. ; qgpath = 0. ; qipath = 0.
   precip_r = 0. ; precip_s = 0. ; precip_g = 0. ; precip_i = 0.
   satrdt(:)= 0. ; supsat(:)= 0. ; tcelci(:)= 0. ; supice(:)= 0.
   ifsat(:) = .false.; ifice(:) =.false.
   ni(:) = 1.e3 ; cldf(:) = 1.
   sumsg(:) = 0. ; vtsg(:) = 0.0
   tstepsnow = 0. ; tstepgraupel = 0.
   pcact(:) = 0. ; nsacw(:) = 0. ; ngacw(:) = 0. ; naacw(:) = 0.
   niacr(:) = 0. ; nsacr(:) = 0. ; ngacr(:) = 0. ; nseml(:) = 0.
   ngeml(:) = 0. ; nracw(:) = 0. ; nccol(:) = 0. ; nrcol(:) = 0.
   ncact(:) = 0. ; nraut(:) = 0. ; ncevp(:) = 0. ; denqn(:) = 0.
   vtn(:) = 0. ; diac(:) = 0. ; diar(:) = 0. ; diai(:) = 0.
   qrcon(:)= 0. ; qrconcr(:) = 0. ; taucon(:) = 0.
   qnpath = 0.
!===============================================================================
!
! sedimentation of qr qs qg qi and nr
!
!===============================================================================
!
! fall out for rain
!
   if(lqr) then
     ktop = ktopqr
     call slope_rain(1,kdim,ktop,qrs(:,i,1),den(:,i),denfac(:,i),t(:,i), &
        ncr(:,i,3), &
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:))
     nstep = 1
     do k = ktop, kts, -1
       vtn(k) = pvtrn*rslopeb(k,1)*denfac(k,i)
       value = max(vtr(k),vtn(k))
       nstep = max(nstep,ceiling(dtcld/delz(k,i)*value))
     enddo
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vtr(k) = (2*vtr(k)+vtr(k+1)+vtr(k-1))/4.
         vtn(k) = (2*vtn(k)+vtn(k+1)+vtn(k-1))/4.
       enddo
     endif
!
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
         if(qrs(k,i,1)>qrmin) then
           denqr(k) = dend(k,i)*qrs(k,i,1)
           denqn(k) = dend(k,i)*ncr(k,i,3)
         else
           denqr(k) = 0.0
           denqn(k) = 0.0
           vtr(k) = 0.0
           vtn(k) = 0.0
         endif
       enddo
!
       call semi_lagrange_sedim(1,kdim,ktop,dend(:,i),denfac(:,i),t(:,i), &
            delz(:,i),vtr(:),denqr(:),qrpath,dtcfl,lat,i)
       call semi_lagrange_sedim(1,kdim,ktop,dend(:,i),denfac(:,i),t(:,i), &
            delz(:,i),vtn(:),denqn(:),qnpath,dtcfl,lat,i)
       do k = ktop, kts, -1
         if(denqr(k)>qmin) then
           qrs(k,i,1) = max(denqr(k)/dend(k,i),0.)
           ncr(k,i,3) = max(denqn(k)/dend(k,i),0.)
           ncr(k,i,3) = max(ncr(k,i,3),ncr5mm*den(k,i)*qrs(k,i,1)) ! 5 mm
           ncr(k,i,3) = min(ncr(k,i,3),ncr82micro*den(k,i)*qrs(k,i,1)) ! 82 micro
         else
           qrs(k,i,1) = 0.0
           ncr(k,i,3) = 0.0
         endif
       enddo
!
       precip_r = qrpath/dtcld + precip_r ! [kgm-2s-1]
!
       call slope_rain(1,kdim,ktop,qrs(:,i,1),den(:,i),denfac(:,i),t(:,i), &
            ncr(:,i,3), &
            rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:))
       do k = ktop, kts, -1
         vtn(k) = pvtrn*rslopeb(k,1)*denfac(k,i)
       enddo
!-------------------------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [hdc 14]
! (v->r or r->v)
       do k = ktop, kts, -1
       if(qrs(k,i,1)>0. .and. zl(k,i)<cbot(i)) then
         supsat(k) = max(q(k,i),qmin)-qsat_mul(k)
         satrdt(k) = supsat(k)/dtcfl
         coeres = rslope(k,1)*sqrt(rslope(k,1)*rslopeb(k,1))
         prevp(k) = (rh_mul(k)-1.)*ncr(k,i,3)*(precr1*rslope(k,1) &
                   +precr2*venfac(k)*coeres)/ab_mul(k)
         if(prevp(k)<=0.) then
           prevp(k) = max(prevp(k),-qrs(k,i,1)/dtcld)
           prevp(k) = max(prevp(k),satrdt(k)*.5)
!----------------------------------------------------------------
! nrevp: evaporation/condensation rate of rain [lh a14]
! (nr->nccn)
           if(prevp(k)==-qrs(k,i,1)/dtcld) then
             ncr(k,i,1) = ncr(k,i,1) + ncr(k,i,3)/float(niter)
             ncr(k,i,3) = ncr(k,i,3) * (1.-1./float(niter))
           endif
         else
           prevp(k) = min(prevp(k),satrdt(k)*.5)
         endif
         value = - xlv(k,i)*prevp(k)
         t(k,i) = t(k,i) - value/cpm(k,i)*dtcfl
         q(k,i) = q(k,i)-prevp(k)*dtcfl
         qrs(k,i,1) = max(qrs(k,i,1) + prevp(k)*dtcfl,0.)
       endif
       enddo
       prevp(:) = 0.
     enddo
   endif
!-------------------------------------------------------------------------------
! fall out for snow and graupel : mass-weighted combined sedimentation (dhl )
!
   if(lqs .or. lqg) then
     ktop = max(ktopqs,ktopqg)
     call slope_snow(1,kdim,ktop,qrs(:,i,2),den(:,i),denfac(:,i),t(:,i), &
          rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:))
     call slope_graupel(1,kdim,ktop,qrs(:,i,3),den(:,i),denfac(:,i),t(:,i), &
          rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:))
!
     do k = ktop, kts, -1
       sumsg(k) = max( (qrs(k,i,2) + qrs(k,i,3)), qmin)
       if(sumsg(k)>qmin) then
         vtsg(k) = (vts(k)*qrs(k,i,2) + vtg(k)*qrs(k,i,3))/sumsg(k)
       else
         vtsg(k) = 0.
       endif
     enddo
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vtsg(k) = (2*vtsg(k)+vtsg(k+1)+vtsg(k-1))/4.
       enddo
     endif
!
     nstep = max(1,ceiling(maxval(dtcld/delz(:,i)*vtsg(:))))
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
         denqs(k) = dend(k,i)*qrs(k,i,2)
         denqg(k) = dend(k,i)*qrs(k,i,3)
       enddo
!
       call semi_lagrange_sedim(1,kdim,ktop,dend(:,i),denfac(:,i),t(:,i), &
            delz(:,i),vtsg(:),denqs(:),qspath,dtcfl,lat,i)
       call semi_lagrange_sedim(1,kdim,ktop,dend(:,i),denfac(:,i),t(:,i), &
            delz(:,i),vtsg(:),denqg(:),qgpath,dtcfl,lat,i)
       do k = kts, ktop
         qrs(k,i,2) = max(denqs(k)/dend(k,i),0.)
         qrs(k,i,3) = max(denqg(k)/dend(k,i),0.)
       enddo
!
       precip_s = qspath/dtcld + precip_s ! [kgm-2s-1]
       precip_g = qgpath/dtcld + precip_g ! [kgm-2s-1]
!
       call slope_snow(1,kdim,ktop,qrs(:,i,2),den(:,i),denfac(:,i),t(:,i), &
            rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:))
       call slope_graupel(1,kdim,ktop,qrs(:,i,3),den(:,i),denfac(:,i),t(:,i), &
            rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:))
!
       sumsg(:) = 0.0
       do k = ktop, kts, -1
         sumsg(k) = max( (qrs(k,i,2) + qrs(k,i,3)), qmin)
         if(sumsg(k)>qmin) then
           vtsg(k) = (vts(k)*qrs(k,i,2) + vtg(k)*qrs(k,i,3)) /sumsg(k)
         else
           vtsg(k) = 0.
         endif
       enddo
     enddo
   endif
!-------------------------------------------------------------------------------
! fall out for cloud ice
!
   if(lqi) then
     ktop = ktopqi
!-------------------------------------------------------------------------------
! ni: ice crystal number concentraiton [hdc 5c]
!
     do k = ktop, kts, -1
       temp = (den(k,i)*max(qci(k,i,2),qcmin))
       temp = sqrt(sqrt(temp*temp*temp))
       ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
       if(qci(k,i,2)<=0.) then
         vti(k) = 0.
       else
         massi = den(k,i)*qci(k,i,2)/ni(k)
         diai(k) = max(min( dicon*sqrt(massi), dimax), qmin)
         vti(k) = 1.49e4*exp(log(diai(k))*(1.31))
       endif
     enddo
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vti(k) = (2*vti(k)+vti(k+1)+vti(k-1))/4.
       enddo
     endif
!
     nstep = max(1,ceiling(maxval(dtcld/delz(:,i)*vti(:))))
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
         denqi(k) = dend(k,i)*qci(k,i,2)
       enddo
!
       call semi_lagrange_sedim(1,kdim,ktop,dend(:,i),denfac(:,i),t(:,i), &
            delz(:,i),vti(:),denqi(:),qipath,dtcfl,lat,i)
       do k = kts, ktop
         qci(k,i,2) = max(denqi(k)/dend(k,i),0.)
       enddo
!
       precip_i = qipath/dtcld + precip_i ! [kgm-2s-1]
!
       do k = ktop, kts, -1
         temp = (den(k,i)*max(qci(k,i,2),qcmin))
         temp = sqrt(sqrt(temp*temp*temp))
         ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
         if(qci(k,i,2)<=0.) then
           vti(k) = 0.
         else
           massi = den(k,i)*qci(k,i,2)/ni(k)
           diai(k) = max(min(dicon * sqrt(massi),dimax), qmin)
           vti(k) = 1.49e4*exp(log(diai(k))*(1.31))
         endif
       enddo
     enddo
   endif
!===============================================================================
!
! precip (den*qrsi*dz/dt) : [kgm-2s-1] ==> rain ( precip/denr*dt*1000 ) : [mm ]
! for wrf unit is mm, whereas it is m for ufs
!
!===============================================================================
   temp = 1.
   precip_sum = precip_r + precip_s + precip_i + precip_g
   precip_ice = precip_s + precip_i
!
   if(precip_sum>0.) then
     value = precip_sum/denr*dtcld*temp
     rain(i) = value + rain(i)
   endif
!
   if(precip_ice>0.) then
     value = precip_ice/denr*dtcld*temp
     tstepsnow = value + tstepsnow
     if ( present (snow)) then
       snow(i) = value + snow(i)
      endif
   endif
!
   if(precip_g>0.) then
     value = precip_g/denr*dtcld*temp
     tstepgraupel = value + tstepgraupel
     if ( present (graupel)) then
       graupel(i) = value + graupel(i)
      endif
   endif
!
!-------------------------------------------------------------------------------
! obtain in-cloud properties
!
   ktop = max(ktopqc,ktopqi)
   call cldf_mps_diag(1,kdim,t(:,i),p(:,i),q(:,i),qci(:,i,1),qci(:,i,2), &
        dxmeter(i),cldf(:),ktop)
!
   do k = ktop, kts, -1
       ! change condensate variables to in-cloud variables
       if(cldf(k)>0.) then
         qci(k,i,1) = qci(k,i,1) / cldf(k)
         qci(k,i,2) = qci(k,i,2) / cldf(k)
         qrs(k,i,1) = qrs(k,i,1) / cldf(k)
         qrs(k,i,2) = qrs(k,i,2) / cldf(k)
         qrs(k,i,3) = qrs(k,i,3) / cldf(k)
       endif
   enddo
!-------------------------------------------------------------------------------
! update the slope parameters for microphysics computation
!
   ktop = ktopqr
   call slope_rain(1,kdim,ktop,qrs(:,i,1),den(:,i),denfac(:,i),t(:,i), &
        ncr(:,i,3), &
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:))
   ktop = ktopqs
   call slope_snow(1,kdim,ktop,qrs(:,i,2),den(:,i),denfac(:,i),t(:,i), &
        rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:))
   ktop = ktopqg
   call slope_graupel(1,kdim,ktop,qrs(:,i,3),den(:,i),denfac(:,i),t(:,i), &
        rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:))
   ktop = ktopqc
   call slope_cloud(1,kdim,ktop,qci(:,i,1),ncr(:,i,2),den(:,i),denfac(:,i), &
        t(:,i),qcmin,rslopec(:),rslopec2(:),rslopec3(:))
!===============================================================================
!
! melting/freezing after sedimentation
!
!===============================================================================
   ktop = ktopmax
   do k = kts, ktop
     tcelci(k) = t(k,i) - t0c
     if(t(k,i)<t0c) ifice(k) = .true.
     n0sfac(k) = max(min(exp(-alpha*tcelci(k)),n0smax/n0s),1.)
   enddo
!-------------------------------------------------------------------------------
! psmlt: melting of snow [hl a33]
! (t>t0: s->r)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,2)>0.) then
         coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
         psmlt(k) = (1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i))) &
               /((t(k,i))+120.)/(den(k,i)) )*(den(k,i))) &
               /xlf0*(t(k,i)-t0c)*pi/2. &
               *n0sfac(k)*(precs1*rslope2(k,2)+precs2 &
               *venfac(k)*coeres)/den(k,i)
         psmlt(k) = max(min(psmlt(k)*dtcld,qrs(k,i,2)),0.)
!-------------------------------------------------------------------------------
! nsmlt: melting of snow [LH A27]
! (T>T0: ->NR)
         if(qrs(k,i,2)>qrmin) then
           sfac = rslope(k,2)*n0s*n0sfac(k)/qrs(k,i,2)
           ncr(k,i,3) = ncr(k,i,3) + sfac*psmlt(k)
         endif
         qrs(k,i,2) = qrs(k,i,2) - psmlt(k)
         qrs(k,i,1) = qrs(k,i,1) + psmlt(k)
         t(k,i) = t(k,i) - xlf0/cpm(k,i)*psmlt(k)*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgmlt: melting of graupel [hl a23]
! (t>t0: g->r)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,3)>0.) then
         coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
         pgmlt(k) = (1.414e3*(1.496e-6*((t(k,i))*sqrt(t(k,i))) &
                    /((t(k,i))+120.)/(den(k,i)) )*(den(k,i))) &
                    /xlf0*(t(k,i)-t0c)*(precg1*rslope2(k,3) &
                    +precg2*venfac(k)*coeres)/den(k,i)
         pgmlt(k) = max(min(pgmlt(k)*dtcld,qrs(k,i,3)),0.)
!-------------------------------------------------------------------------------
! ngmlt: melting of graupel [lh a28]
! (t>t0: ->nr)
       if(qrs(k,i,3)>qrmin) then
         gfac = rslope(k,3)*n0g/qrs(k,i,3)
         ncr(k,i,3) = ncr(k,i,3) + gfac*pgmlt(k)
       endif
         qrs(k,i,3) = qrs(k,i,3) - pgmlt(k)
         qrs(k,i,1) = qrs(k,i,1) + pgmlt(k)
         t(k,i) = t(k,i) - xlf0/cpm(k,i)*pgmlt(k)*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [hl a47]
! (t>t0: i->c)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qci(k,i,2)>0.) then
         qci(k,i,1) = qci(k,i,1) + qci(k,i,2)
         t(k,i) = t(k,i) - xlf0/cpm(k,i)*qci(k,i,2)*cldf(k)
         qci(k,i,2) = 0.
!---------------------------------------------------------------
! nimlt: instantaneous melting of cloud ice [lh a18]
! (t>t0: ->nc)
         ncr(k,i,2) = ncr(k,i,2) + ni(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [hl a45]
! (t<-40c: c->i)
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
       if(tcelci(k)<-40. .and. qci(k,i,1)>0.) then
         qci(k,i,2) = qci(k,i,2) + qci(k,i,1)
         t(k,i) = t(k,i) + xlf(k,i)/cpm(k,i)*qci(k,i,1)*cldf(k)
         qci(k,i,1) = 0.
!---------------------------------------------------------------
! nihmf: homogeneous freezing of cloud water below -40c [lh a17]
! (t<-40c: nc-> )
         if(ncr(k,i,2)>0.) ncr(k,i,2) = 0.
       endif
!-------------------------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [hl a44]
! (t0>t>-40c: c->i)
       if(ifice(k) .and. qci(k,i,1)>qcmin) then
         value = max(tcelci(k),-70.)
         pfrzdtc = min(pisq*pfrz1*(exp(-pfrz2*value)-1.)*denr/den(k,i) &
                   *ncr(k,i,2)*rslopec3(k)*rslopec3(k)/18.*dtcld &
                   ,qci(k,i,1))
!---------------------------------------------------------------
! nihtf: heterogeneous freezing of cloud water [LH A16]
! (T0>T>-40C: NC->)
         if(ncr(k,i,2)>ncmin) then
           nfrzdtc = min(pi*pfrz1*(exp(-pfrz2*value)-1.)*ncr(k,i,2) &
                     *rslopec3(k)/6.*dtcld,ncr(k,i,2))
           ncr(k,i,2) = max(ncr(k,i,2) - nfrzdtc,0.)
         endif
         qci(k,i,1) = qci(k,i,1) - pfrzdtc
         qci(k,i,2) = qci(k,i,2) + pfrzdtc
         t(k,i) = t(k,i) + xlf(k,i)/cpm(k,i)*pfrzdtc*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgfrz: freezing of rain water [hl a20]
! (t<t0, r->g)
   if(lqr) then
     ktop = ktopqr
     do k = ktop, kts, -1
       if(ifice(k) .and. qrs(k,i,1)>0.) then
         value = max(tcelci(k), -70.)
         pfrzdtr = min(140.*pisq*pfrz1*ncr(k,i,3)*denr/den(k,i) &
                  *(exp(-pfrz2*value)-1.)*rslope3(k,1)*rslope3(k,1) &
                  *dtcld,qrs(k,i,1))
!---------------------------------------------------------------
! ngfrz: freezing of rain water [lh a26]
! (T<T0, NR-> )
         if(ncr(k,i,3)>nrmin) then
           nfrzdtr = min(4.*pi*pfrz1*ncr(k,i,3)*(exp(-pfrz2*value)-1.) &
                    *rslope3(k,1)*dtcld, ncr(k,i,3))
           ncr(k,i,3) = max(ncr(k,i,3) - nfrzdtr, 0.0)
         endif
         qrs(k,i,1) = qrs(k,i,1) - pfrzdtr
         qrs(k,i,3) = qrs(k,i,3) + pfrzdtr
         t(k,i) = t(k,i) + xlf(k,i)/cpm(k,i)*pfrzdtr*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! re-define cloud top for numerical efficiencies
!
   ktop = ktopini
   do k = kts, ktop
     if(t(k,i)<ttp) then
       xal = xai
       xbl = xbi
     else
       xal = xa
       xbl = xb
     endif
     tr=ttp/t(k,i)
     logtr=log(tr)
     qsat_mul(k) = psat*exp(logtr*(xa)+xb*(1.-tr))
     qsat_mul(k) = min(qsat_mul(k),0.99*p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k) = max(q(k,i) / qsat_mul(k),qmin)
     qsat_ice(k) = psat*exp(logtr*(xal)+xbl*(1.-tr))
     qsat_ice(k) = min(qsat_ice(k),0.99*p(k,i))
     qsat_ice(k) = ep2 * qsat_ice(k) / (p(k,i) - qsat_ice(k))
     qsat_ice(k) = max(qsat_ice(k),qmin)
     rh_ice(k) = max(q(k,i) / qsat_ice(k),qmin)
   enddo
!
   lqc = .false.
   lqi = .false.
   lqr = .false.
   lqs = .false.
   lqg = .false.
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,qci(:,i,2),zero_0,ktopqi)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,2),zero_0,ktopqs)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,3),zero_0,ktopqg)
   call find_cloud_top(1,kdim,ktopini,rh_ice(:), one_1,ktoprh)
   if(ktopqc>0.0) lqc = .true.
   if(ktopqi>0.0) lqi = .true.
   if(ktopqr>0.0) lqr = .true.
   if(ktopqs>0.0) lqs = .true.
   if(ktopqg>0.0) lqg = .true.
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktoprh)
!-------------------------------------------------------------------------------
! restore grid-mean variables
!
   ktop = max(ktopqc,ktopqi)
   do k = ktop, kts, -1
       ! change in-cloud condensate variables to grid-mean variables
     if(cldf(k)>0.) then
       qci(k,i,1) = qci(k,i,1) * cldf(k)
       qci(k,i,2) = qci(k,i,2) * cldf(k)
       qrs(k,i,1) = qrs(k,i,1) * cldf(k)
       qrs(k,i,2) = qrs(k,i,2) * cldf(k)
       qrs(k,i,3) = qrs(k,i,3) * cldf(k)
! ncr(k,i,2) = ncr(k,i,2) * cldf(k)
! ncr(k,i,3) = ncr(k,i,3) * cldf(k)
     endif
   enddo
!-------------------------------------------------------------------------------
! obtain in-cloud properties
!
   call cldf_mps_diag(1,kdim,t(:,i),p(:,i),q(:,i),qci(:,i,1),qci(:,i,2), &
        dxmeter(i),cldf(:),ktop)
   do k = ktop, kts, -1
       ! change condensate variables to in-cloud variables
     if(cldf(k)>0.) then
       qci(k,i,1) = qci(k,i,1) / cldf(k)
       qci(k,i,2) = qci(k,i,2) / cldf(k)
       qrs(k,i,1) = qrs(k,i,1) / cldf(k)
       qrs(k,i,2) = qrs(k,i,2) / cldf(k)
       qrs(k,i,3) = qrs(k,i,3) / cldf(k)
     endif
   enddo
!-------------------------------------------------------------------------------
! update the slope parameters and microphysical parameters
!
   ktop = ktopqr
   call slope_rain(1,kdim,ktop,qrs(:,i,1),den(:,i),denfac(:,i),t(:,i), &
        ncr(:,i,3), &
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:))
   ktop = ktopqs
   call slope_snow(1,kdim,ktop,qrs(:,i,2),den(:,i),denfac(:,i),t(:,i), &
        rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:))
   ktop = ktopqg
   call slope_graupel(1,kdim,ktop,qrs(:,i,3),den(:,i),denfac(:,i),t(:,i), &
        rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:))
   ktop = ktopqi
   do k = ktop, kts, -1
     temp = (den(k,i)*max(qci(k,i,2),qcmin))
     temp = sqrt(sqrt(temp*temp*temp))
     ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
     massi = den(k,i)*qci(k,i,2)/ni(k)
     diai(k) = max(min(dicon * sqrt(massi),dimax), qmin)
     vti(k) = 1.49e4*exp(log(diai(k))*(1.31))
   enddo
!
   ktop = max(ktopqs,ktopqg)
   do k = ktop, kts, -1
     sumsg(k) = max( (qrs(k,i,2)+qrs(k,i,3)), qmin)
     if(sumsg(k)>qmin) then
       vtsg(k) = (vts(k)*qrs(k,i,2)+vtg(k)*qrs(k,i,3))/sumsg(k)
     else
       vtsg(k) = 0.
     endif
   enddo
!-------------------------------------------------------------------------------
! compute the mean-volume drop diameter for raindrop distribution [lh a10]
!
   ktop = ktopqr
   do k = ktop, kts, -1
     diar(k) = rslope(k,1)*drcoeff
   enddo
!--------------------------------------------------------------------
! compute the mean-volume drop diameter for cloud-droplet distribution [lh a7]
!
   ktop = ktopqc
   call slope_cloud(1,kdim,ktop,qci(:,i,1),ncr(:,i,2),den(:,i),denfac(:,i), &
        t(:,i),qcmin,rslopec(:),rslopec2(:),rslopec3(:))
   do k = ktop, kts, -1
     diac(k) = rslopec(k)
   enddo
!===============================================================================
!
! warm rain processes
!
!===============================================================================
   ktop = max(ktopqc,ktopqr)
   do k = ktop, kts, -1
     qrcon(k) = 2.7e-2*den(k,i)*qci(k,i,1)*(1.e20/16.*rslopec2(k) &
               *rslopec2(k)-0.4)
     qrconcr(k) = max(1.2*qrcon(k), qrmin)
   enddo
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [lh 9]
! (qc->qr)
       if(diac(k)>di15) then
         taucon(k) = 3.7/den(k,i)/qci(k,i,1)/(0.5e6*rslopec(k)-7.5)
         taucon(k) = max(taucon(k), qmin)
         praut(k) = qrcon(k)/(taucon(k)*den(k,i))
         praut(k) = min(max(praut(k),0.),qci(k,i,1)*rdtcld)
!---------------------------------------------------------------
! nraut: auto conversion rate from cloud to rain [lh a6] [cp 18 & 19]
! (nc->nr)
         nraut(k) = 3.5e9*den(k,i)*praut(k)
         if(qrs(k,i,1)>qrconcr(k)) then
           nraut(k) = ncr(k,i,3)/qrs(k,i,1)*praut(k)
           nraut(k) = min(nraut(k),ncr(k,i,2)*rdtcld)
         endif
       endif
     enddo
   endif
!
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
!----------------------------------------------------------------
! nccol: self collection of cloud water [lh a8] [cp 24 & 25]
! (nc->)
       if(diac(k)>=di100) then
         nccol(k) = ncrk1*ncr(k,i,2)*ncr(k,i,2)*rslopec3(k)
       else
         nccol(k) = 2.*ncrk2*ncr(k,i,2)*ncr(k,i,2)*rslopec3(k)*rslopec3(k)
       endif
     enddo
   endif
!
   if(lqr) then
     ktop = ktopqr
     do k = ktop, kts, -1
       supsat(k) = max(q(k,i),qmin)-qsat_mul(k)
       satrdt(k) = supsat(k)*rdtcld
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [lh 10] [cp 22 & 23]
! (qc->qr)
! nracw: accretion of cloud water by rain [lh a9]
! (nc->)
       if(qrs(k,i,1)>=qrconcr(k)) then
         if(diar(k)>=di100) then
           nracw(k) = min(ncrk1*ncr(k,i,2)*ncr(k,i,3)*(rslopec3(k) &
                      + 24.*rslope3(k,1)),ncr(k,i,2)*rdtcld)
           pracw(k) = min(pi/6.*(denr/den(k,i))*ncrk1*ncr(k,i,2) &
                      *ncr(k,i,3)*rslopec3(k)*(2.*rslopec3(k) &
                      + 24.*rslope3(k,1)),qci(k,i,1)*rdtcld)
         else
           nracw(k) = min(ncrk2*ncr(k,i,2)*ncr(k,i,3)*(2.*rslopec3(k) &
                      *rslopec3(k)+5040.*rslope3(k,1) &
                      *rslope3(k,1)),ncr(k,i,2)*rdtcld)
           pracw(k) = min(pi/6.*(denr/den(k,i))*ncrk2*ncr(k,i,2) &
                      *ncr(k,i,3)*rslopec3(k)*(6.*rslopec3(k) &
                      *rslopec3(k)+5040.*rslope3(k,1)*rslope3(k,1)) &
                      ,qci(k,i,1)*rdtcld)
         endif
       endif
!----------------------------------------------------------------
! nrcol: self collection of rain-drops and break-up [lh a21] [cp 24 & 25]
! (nr->)
       if(qrs(k,i,1)>=qrconcr(k)) then
         if(diar(k)<di100) then
           nrcol(k) = 5040.*ncrk2*ncr(k,i,3)*ncr(k,i,3)*rslope3(k,1) &
                     *rslope3(k,1)
         elseif(diar(k)>=di100 .and. diar(k)<di600) then
           nrcol(k) = 24.*ncrk1*ncr(k,i,3)*ncr(k,i,3)*rslope3(k,1)
         elseif(diar(k)>=di600 .and. diar(k)<di2000) then
           value = -2.5e3*(diar(k)-di600)
           nrcol(k) = 24.*exp(value)*ncrk1*ncr(k,i,3)*ncr(k,i,3) &
                      *rslope3(k,1)
         else
           nrcol(k) = 0.
         endif
       endif
!-------------------------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [hdc 14]
! (v->r or r->v)
       if(qrs(k,i,1)>0.) then
         coeres = rslope(k,1)*sqrt(rslope(k,1)*rslopeb(k,1))
         prevp(k) = (rh_mul(k)-1.)*ncr(k,i,3)*(precr1*rslope(k,1) &
                   +precr2*venfac(k)*coeres)/ab_mul(k)
         if(prevp(k)<=0.) then
           prevp(k) = max(prevp(k),-qrs(k,i,1)*rdtcld)
           prevp(k) = max(prevp(k),satrdt(k)*.5)
!----------------------------------------------------------------
! nrevp: evaporation/condensation rate of rain [lh a14]
! (nr->nccn)
           if(prevp(k)==-qrs(k,i,1)*rdtcld) then
             ncr(k,i,1) = ncr(k,i,1) + ncr(k,i,3)
             ncr(k,i,3) = 0.
           endif
         else
           prevp(k) = min(prevp(k),satrdt(k)*.5)
         endif
       endif
! prevp(k) = 0.
     enddo
   endif
!===============================================================================
!
! cold rain processes
!
!===============================================================================
   ktop = ktopmax
   ifice(:) = .false.
   do k = ktop, kts, -1
     tcelci(k) = t(k,i) - t0c
     if(t(k,i)<t0c) ifice(k) = .true.
     n0sfac(k) = max(min(exp(-alpha*tcelci(k)),n0smax/n0s),1.)
     supsat(k) = max(q(k,i),qmin)-qsat_ice(k)
     satrdt(k) = supsat(k)*rdtcld
   enddo
!
   if(lqr .and. lqi) then
     ktop = min(ktopqr,ktopqi)
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,1)>qrmin .and. qci(k,i,2)>qcmin) then
!-------------------------------------------------------------
! praci: accretion of cloud ice by rain [hl a15] [lfo 25]
! (t<t0: qi->qr)
           acrfac = 6.*rslope2(k,1)+4.*diai(k)*rslope(k,1) + diai(k)**2
           praci(k) = pi*qci(k,i,2)*ncr(k,i,3)*abs(vtr(k)-vti(k))*acrfac/4.
             ! reduce collection efficiency (suggested by B. Wilt)
           praci(k) = praci(k)*min(max(0.0,qrs(k,i,1)/qci(k,i,2)),1.)**2
           praci(k) = min(praci(k),qci(k,i,2)*rdtcld)
!-------------------------------------------------------------
! piacr: accretion of rain by cloud ice [hl A19] [lfo 26]
! (t<t0: qr->qs OR qr->qg)
           piacr(k) = pisq*avtr*ncr(k,i,3)*denr*ni(k)*denfac(k,i) &
                       *g7pbr*rslope3(k,1)*rslope2(k,1)*rslopeb(k,1) &
                       /24./den(k,i)
             ! reduce collection efficiency (suggested by B. Wilt)
           piacr(k) = piacr(k)*min(max(0.0,qci(k,i,2)/qrs(k,i,1)),1.)**2
           piacr(k) = min(piacr(k),qrs(k,i,1)*rdtcld)
!-------------------------------------------------------------
! niacr: accretion of rain by cloud ice [lh a25]
! (t<t0: nr->)
           if(ncr(k,i,3)>nrmin) then
             niacr(k) = pi*avtr*ncr(k,i,3)*ni(k)*denfac(k,i)*g4pbr &
                         *rslope2(k,1)*rslopeb(k,1)/4.
             ! reduce collection efficiency (suggested by B. Wilt)
             niacr(k) = niacr(k)*min(max(0.0,qci(k,i,2)/qrs(k,i,1)),1.)**2
             niacr(k) = min(niacr(k),ncr(k,i,3)*rdtcld)
           endif
         endif
       endif
     enddo
   endif
!
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
!
       if(ifice(k)) then
         if(qrs(k,i,2)>qrmin .and. qci(k,i,2)>qcmin) then
           eacrs = exp(0.09*(tcelci(k)))
!-------------------------------------------------------------------------------
! psaci: accretion of cloud ice by rain [hdc 10]
! (t<t0: i->s)
           acrfac = 2.*rslope3(k,2) + 2.*diai(k)*rslope2(k,2) &
                   + diai(k)**2*rslope(k,2)
           psaci(k) = pi*qci(k,i,2)*eacrs*n0s*n0sfac(k) &
                        *abs(vtsg(k)-vti(k))*acrfac/4.
           psaci(k) = psaci(k)*min(max(0.0,qrs(k,i,2)/qci(k,i,2)),1.)**2
           psaci(k) = min(psaci(k),qci(k,i,2)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! psacw: accretion of cloud water by snow [hl a7]
! (t<t0: c->s, and t>=t0: c->r)
       if(qrs(k,i,2)>qrmin .and. qci(k,i,1)>qcmin) then
         psacw(k) = min(pacrc*n0sfac(k)*rslope3(k,2)*rslopeb(k,2) &
                  *qci(k,i,1)*denfac(k,i),qci(k,i,1)*rdtcld)
             ! reduce collection efficiency (suggested by B. Wilt)
         psacw(k) = psacw(k)*min(max(0.0,qrs(k,i,2)/qci(k,i,1)),1.)**2
         psacw(k) = min(psacw(k),qci(k,i,1)*rdtcld)
       endif
!-------------------------------------------------------------
! nsacw: accretion of cloud water by snow [lh A12]
! (NC ->)
       if(qrs(k,i,2)>qrmin .and. ncr(k,i,2)>ncmin) then
         nsacw(k) = min(pacrc*n0sfac(k)*rslope3(k,2)*rslopeb(k,2) &
             ! reduce collection efficiency (suggested by B. Wilt)
                     *min(max(0.0,qrs(k,i,2)/qci(k,i,1)),1.)**2 &
                     *ncr(k,i,2)*denfac(k,i),ncr(k,i,2)*rdtcld)
       endif
     enddo
   endif
!
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(ifice(k)) then
!-------------------------------------------------------------------------------
! pgaci: accretion of cloud ice by graupel [hl a17]
! (t<t0: i->g)
         if(qrs(k,i,3)>qrmin .and. qci(k,i,2)>qcmin) then
           egi = exp(0.09*(tcelci(k)))
           acrfac = 2.*rslope3(k,3) + 2.*diai(k)*rslope2(k,3) &
                    + diai(k)**2*rslope(k,3)
           pgaci(k) = pi*egi*qci(k,i,2)*n0g*abs(vtsg(k)-vti(k))*acrfac/4.
           pgaci(k) = min(pgaci(k),qci(k,i,2)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! pgacw: accretion of cloud water by graupel [hl a6]
! (t<t0: c->g, and t>=t0: c->r)
       if(qrs(k,i,3)>qrmin .and. qci(k,i,1)>qcmin) then
          pgacw(k) = min(pacrg*rslope3(k,3)*rslopeb(k,3) &
             ! reduce collection efficiency (suggested by B. Wilt)
                    *min(max(0.0,qrs(k,i,3)/qci(k,i,1)),1.)**2 &
                    *qci(k,i,1)*denfac(k,i),qci(k,i,1)*rdtcld)
       endif
!-------------------------------------------------------------
! ngacw: accretion of cloud water by graupel [lh a13]
! (nc->
       if(qrs(k,i,3)>qrmin .and. ncr(k,i,2)>ncmin) then
         ngacw(k) = min(pacrg*rslope3(k,3)*rslopeb(k,3)*ncr(k,i,2) &
             ! reduce collection efficiency (suggested by B. Wilt)
                     *min(max(0.0,qrs(k,i,3)/qci(k,i,1)),1.)**2 &
                     *denfac(k,i),ncr(k,i,2)*rdtcld)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! paacw: accretion of cloud water by averaged snow/graupel
! (t<t0: c->g or s, and t>=t0: c->r)
   ktop = max(ktopqs,ktopqg)
   do k = ktop, kts, -1
     if(sumsg(k)>qmin) then
       paacw(k) = (qrs(k,i,2)*psacw(k) + qrs(k,i,3)*pgacw(k))/sumsg(k)
!-------------------------------------------------------------
! naacw: Accretion of cloud water by averaged snow/graupel
! (nc->)
       naacw(k) = (qrs(k,i,2)*nsacw(k) + qrs(k,i,3)*ngacw(k))/sumsg(k)
     endif
   enddo
!
   if(lqr) then
     ktop = ktopqr
     do k = kts, ktop
!-------------------------------------------------------------------------------
! pracs: accretion of snow by rain [hl a11]
! (t<t0: s->g)
       if(qrs(k,i,2)>qrmin .and. qrs(k,i,1)>qrmin) then
         if(ifice(k)) then
           acrfac = 5.*rslope3(k,2)*rslope3(k,2) &
                  + 4.*rslope3(k,2)*rslope2(k,2)*rslope(k,1) &
                  + 1.5*rslope2(k,2)*rslope2(k,2)*rslope2(k,1)
           pracs(k) = pisq*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtr(k)-vtsg(k)) &
                      *(dens/den(k,i))*acrfac
             ! reduce collection efficiency (suggested by B. Wilt)
           pracs(k) = pracs(k)*min(max(0.0,qrs(k,i,1)/qrs(k,i,2)),1.)**2
           pracs(k) = min(pracs(k),qrs(k,i,2)*rdtcld)
         endif
!-------------------------------------------------------------
! psacr: accretion of rain by snow [hl a10] [lfo 28]
! (t<t0:qr->qs or qr->qg) (t>=t0: enhanced melting of snow)
         acrfac = 30.*rslope3(k,1)*rslope2(k,1)*rslope(k,2) &
                 +10.*rslope2(k,1)*rslope2(k,1)*rslope2(k,2) &
                 + 2.*rslope3(k,1)*rslope3(k,2 )
         psacr(k) = pisq*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtsg(k)-vtr(k)) &
                        *(denr/den(k,i))*acrfac
             ! reduce collection efficiency (suggested by B. Wilt)
         psacr(k) = psacr(k)*min(max(0.0,qrs(k,i,2)/qrs(k,i,1)),1.)**2
         psacr(k) = min(psacr(k),qrs(k,i,1)*rdtcld)
       endif
       if(qrs(k,i,2)>qrmin .and. ncr(k,i,3)>nrmin) then
!-------------------------------------------------------------
! nsacr: accretion of rain by snow [LH A23]
! (t<t0:nr->)
           acrfac = 1.5*rslope2(k,1)*rslope(k,2) &
                  + 1.0*rslope(k,1)*rslope2(k,2) + .5*rslope3(k,2)
           nsacr(k) = pi*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtsg(k)-vtr(k)) &
                        *acrfac
             ! reduce collection efficiency (suggested by B. Wilt)
           nsacr(k) = nsacr(k)*min(max(0.0,qrs(k,i,2)/qrs(k,i,1)),1.)**2
           nsacr(k) = min(nsacr(k),ncr(k,i,3)*rdtcld)
       endif
!-------------------------------------------------------------------------------
! pgacr: accretion of rain by graupel [hl a12]
! (t<t0: r->g) (t>=t0: enhanced melting of graupel)
       if(qrs(k,i,3)>qrmin .and. qrs(k,i,1)>qrmin) then
!-------------------------------------------------------------
! pgacr: accretion of rain by graupel [HL A12] [LFO 42]
! (t<t0: qr->qg) (t>=t0: enhanced melting of graupel)
         acrfac = 30.*rslope3(k,1)*rslope2(k,1)*rslope(k,3) &
                 + 10.*rslope2(k,1)*rslope2(k,1)*rslope2(k,3) &
                 + 2.*rslope3(k,1)*rslope3(k,3)
         pgacr(k) = pisq*ncr(k,i,3)*n0g*abs(vtsg(k)-vtr(k))*(denr/den(k,i)) &
                       *acrfac
             ! reduce collection efficiency (suggested by B. Wilt)
         pgacr(k) = pgacr(k)*min(max(0.0,qrs(k,i,3)/qrs(k,i,1)),1.)**2
         pgacr(k) = min(pgacr(k),qrs(k,i,1)*rdtcld)
       endif
!-------------------------------------------------------------
! ngacr: Accretion of rain by graupel [LH A24]
! (t<t0: nr->)
       if(qrs(k,i,3)>qrmin .and. ncr(k,i,3)>nrmin) then
         acrfac = 1.5*rslope2(k,1)*rslope(k,3) &
                 +1.0*rslope(k,1)*rslope2(k,3) + .5*rslope3(k,3)
         ngacr(k) = pi*ncr(k,i,3)*n0g*abs(vtsg(k)-vtr(k))*acrfac
           ! reduce collection efficiency (suggested by B. Wilt)
         ngacr(k) = ngacr(k)*min(max(0.0,qrs(k,i,3)/qrs(k,i,1)),1.)**2
         ngacr(k) = min(ngacr(k),ncr(k,i,3)*rdtcld)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pseml: enhanced melting of snow by accretion of water [hl a34]
! (t>=t0: s->r)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,2)>0.) then
         pseml(k) = min(max(-cliq*tcelci(k)*(paacw(k)+psacr(k))/xlf0, &
                    -qrs(k,i,2)*rdtcld),0.)
!--------------------------------------------------------------
! nseml: enhanced melting of snow by accretion of water [Lh a29]
! (t>=t0: ->nr)
         if (qrs(k,i,2)>qrmin) then
           sfac = rslope(k,2)*n0s*n0sfac(k)/qrs(k,i,2)
           nseml(k) = -sfac*pseml(k)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgeml: enhanced melting of graupel by accretion of water [hl a24]
! (t>=t0: g->r)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,3)>0.) then
         pgeml(k) = min(max(-cliq*tcelci(k)*(paacw(k)+pgacr(k))/xlf0, &
                   -qrs(k,i,3)*rdtcld),0.)
!--------------------------------------------------------------
! ngeml: enhanced melting of graupel by accretion of water [lh a30]
! (t>=t0: -> nr)
         if (qrs(k,i,3)>qrmin) then
           gfac = rslope(k,3)*n0g/qrs(k,i,3)
           ngeml(k) = -gfac*pgeml(k)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pidep: deposition/sublimation rate of ice [hdc 9]
! (t<t0: v->i or i->v)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qci(k,i,2)>0 .and. (.not.ifsat(k))) then
           pidep(k) = 4.*diai(k)*ni(k)*(rh_ice(k)-1.)/ab_ice(k)
           supice(k) = satrdt(k) - prevp(k)
           if(pidep(k)<0.) then
             pidep(k) = max(max(pidep(k),satrdt(k)*.5),supice(k))
             pidep(k) = max(pidep(k),-qci(k,i,2)*rdtcld)
           else
             pidep(k) = min(min(pidep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k))>=abs(satrdt(k))) ifsat(k) = .true.
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [hdc 14]
! (v->s or s->v)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,2)>0. .and. (.not.ifsat(k))) then
           coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
           psdep(k) = (rh_ice(k)-1.)*n0sfac(k)*(precs1*rslope2(k,2) &
                      + precs2*venfac(k)*coeres)/ab_ice(k)
           supice(k) = satrdt(k)-prevp(k)-pidep(k)
           if(psdep(k)<0.) then
             psdep(k) = max(psdep(k),-qrs(k,i,2)*rdtcld)
             psdep(k) = max(max(psdep(k),satrdt(k)*.5),supice(k))
           else
             psdep(k) = min(min(psdep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k)+psdep(k))>=abs(satrdt(k))) &
             ifsat(k) = .true.
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [hl a21]
! (t<t0: v->g or g->v)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,3)>0. .and. (.not.ifsat(k))) then
           coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
           pgdep(k) = (rh_ice(k)-1.)*(precg1*rslope2(k,3) &
                     + precg2*venfac(k)*coeres)/ab_ice(k)
           supice(k) = satrdt(k)-prevp(k)-pidep(k)-psdep(k)
           if(pgdep(k)<0.) then
             pgdep(k) = max(pgdep(k),-qrs(k,i,3)*rdtcld)
             pgdep(k) = max(max(pgdep(k),satrdt(k)*.5),supice(k))
           else
             pgdep(k) = min(min(pgdep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k)+psdep(k)+pgdep(k))>=abs(satrdt(k))) &
             ifsat(k) = .true.
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [hl a50] [hdc 7-8]
! (t<t0: v->i)
   ktop = ktoprh
   do k = ktop, kts, -1
     if(ifice(k)) then
       if(supsat(k)>0 .and. (.not.ifsat(k))) then
         supice(k) = satrdt(k)-prevp(k)-pidep(k)-psdep(k)-pgdep(k)
         ni0 = 117.*dust0*exp(-0.125*tcelci(k))
         roqi0 = 4.92e-11*exp(log(ni0)*(1.33))
         pigen(k) = max(0.,(roqi0/den(k,i)-max(qci(k,i,2),0.))*rdtcld)
         pigen(k) = min(min(pigen(k),satrdt(k)),supice(k))
       endif
     endif
   enddo
!-------------------------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [hdc 12]
! (t<t0: i->s)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qci(k,i,2)>0.) then
           qimax = roqimax/den(k,i)
           psaut(k) = max(0.,(qci(k,i,2)-qimax)*rdtcld)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [hl a4]
! (t<t0: s->g)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,2)>0.) then
           alpha2 = 1.e-3*exp(0.09*(tcelci(k)))
           pgaut(k) = min(max(0.,alpha2*(qrs(k,i,2)-qs0)),qrs(k,i,2)*rdtcld)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! psevp: evaporation of melting snow [hl a35]
! (t>t0: s->v)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(.not.ifice(k)) then
         if(qrs(k,i,2)>0. .and. rh_mul(k)<1.) then
              coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
              psevp(k) = (rh_mul(k)-1.)*n0sfac(k)*(precs1*rslope2(k,2) &
                       + precs2*venfac(k)*coeres)/ab_mul(k)
           psevp(k) = min(max(psevp(k),-qrs(k,i,2)*rdtcld),0.)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgevp: evaporation of melting graupel [hl a25]
! (t>=t0: g->v)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(.not.ifice(k)) then
         if(qrs(k,i,3)>0. .and. rh_mul(k)<1.) then
           coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
           pgevp(k) = (rh_mul(k)-1.)*(precg1*rslope2(k,3) &
                     + precg2*venfac(k)*coeres)/ab_mul(k)
           pgevp(k) = min(max(pgevp(k),-qrs(k,i,3)*rdtcld),0.)
         endif
       endif
     enddo
   endif
!===============================================================================
!
! convert the in-cloud source/sink terms to the grd-mean state
!
!===============================================================================
   ktop = max(ktopqc,ktopqi)
   do k = ktop, kts, -1
       if(cldf(k)>0.) then
           ! change in-cloud condensate variables to grid-mean variables
         qci(k,i,1) = qci(k,i,1) * cldf(k)
         qci(k,i,2) = qci(k,i,2) * cldf(k)
         qrs(k,i,1) = qrs(k,i,1) * cldf(k)
         qrs(k,i,2) = qrs(k,i,2) * cldf(k)
         qrs(k,i,3) = qrs(k,i,3) * cldf(k)
           ! get grid-mean rates
         praut(k) = praut(k) * cldf(k)
         pracw(k) = pracw(k) * cldf(k)
         prevp(k) = prevp(k) * cldf(k)
         praci(k) = praci(k) * cldf(k)
         piacr(k) = piacr(k) * cldf(k)
         psaci(k) = psaci(k) * cldf(k)
         pgaci(k) = pgaci(k) * cldf(k)
         psacw(k) = psacw(k) * cldf(k)
         pgacw(k) = pgacw(k) * cldf(k)
         paacw(k) = paacw(k) * cldf(k)
         pracs(k) = pracs(k) * cldf(k)
         pgacr(k) = pgacr(k) * cldf(k)
         pgacs(k) = pgacs(k) * cldf(k)
         pseml(k) = pseml(k) * cldf(k)
         pgeml(k) = pgeml(k) * cldf(k)
         pidep(k) = pidep(k) * cldf(k)
         psdep(k) = psdep(k) * cldf(k)
         pgdep(k) = pgdep(k) * cldf(k)
         pigen(k) = pigen(k) * cldf(k)
         psaut(k) = psaut(k) * cldf(k)
         pgaut(k) = pgaut(k) * cldf(k)
         psevp(k) = psevp(k) * cldf(k)
         pgevp(k) = pgevp(k) * cldf(k)
       endif
   enddo
!===============================================================================
!
! check mass conservation of source/sink terms and feedback to the large scale
!
!===============================================================================
   ktop = ktopmax
   do k = ktop, kts, -1
!
     delta2 = 0.
     delta3 = 0.
     if(qrs(k,i,1)<1.e-4 .and. qrs(k,i,2)<1.e-4) delta2 = 1.
     if(qrs(k,i,1)<1.e-4) delta3 = 1.
!
     if(ifice(k)) then
!
! cloud water
!
       value = max(qmin,qci(k,i,1))
       source = (praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
       if (source>value) then
         factor = value/source
         praut(k) = praut(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
       endif
!
! cloud ice
!
       value = max(qmin,qci(k,i,2))
       source = (psaut(k)-pigen(k)-pidep(k)+praci(k)+psaci(k)+pgaci(k))*dtcld
       if (source>value) then
         factor = value/source
         psaut(k) = psaut(k) * factor
         pigen(k) = pigen(k) * factor
         pidep(k) = pidep(k) * factor
         praci(k) = praci(k) * factor
         psaci(k) = psaci(k) * factor
         pgaci(k) = pgaci(k) * factor
       endif
!
! rain
!
       value = max(qmin,qrs(k,i,1))
       source = (-praut(k)-prevp(k)-pracw(k)+piacr(k)+psacr(k)+pgacr(k))*dtcld
       if (source>value) then
         factor = value/source
         praut(k) = praut(k) * factor
         prevp(k) = prevp(k) * factor
         pracw(k) = pracw(k) * factor
         piacr(k) = piacr(k) * factor
         psacr(k) = psacr(k) * factor
         pgacr(k) = pgacr(k) * factor
       endif
!
! snow
!
       value = max(qmin,qrs(k,i,2))
       source = - (psdep(k)+psaut(k)-pgaut(k)+paacw(k)+piacr(k)*delta3 &
                + praci(k)*delta3 - pracs(k)*(1.-delta2) &
                + psacr(k)*delta2 + psaci(k)-pgacs(k) )*dtcld
       if (source>value) then
         factor = value/source
         psdep(k) = psdep(k) * factor
         psaut(k) = psaut(k) * factor
         pgaut(k) = pgaut(k) * factor
         paacw(k) = paacw(k) * factor
         piacr(k) = piacr(k) * factor
         praci(k) = praci(k) * factor
         psaci(k) = psaci(k) * factor
         pracs(k) = praCs(k) * factor
         psacr(k) = psacr(k) * factor
         pgacs(k) = pgacs(k) * factor
       endif
!
! graupel
!
       value = max(qmin,qrs(k,i,3))
       source = - (pgdep(k)+pgaut(k) &
                + piacr(k)*(1.-delta3) + praci(k)*(1.-delta3) &
                + psacr(k)*(1.-delta2) + pracs(k)*(1.-delta2) &
                + pgaci(k)+paacw(k)+pgacr(k)+pgacs(k))*dtcld
       if (source>value) then
         factor = value/source
         pgdep(k) = pgdep(k) * factor
         pgaut(k) = pgaut(k) * factor
         piacr(k) = piacr(k) * factor
         praci(k) = praci(k) * factor
         psacr(k) = psacr(k) * factor
         pracs(k) = pracs(k) * factor
         paacw(k) = paacw(k) * factor
         pgaci(k) = pgaci(k) * factor
         pgacr(k) = pgacr(k) * factor
         pgacs(k) = pgacs(k) * factor
       endif
!
! cloud
!
       value = max(qmin,ncr(k,i,2))
       source = (nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
       if (source>value) then
         factor = value/source
         nraut(k) = nraut(k) * factor
         nccol(k) = nccol(k) * factor
         nracw(k) = nracw(k) * factor
         naacw(k) = naacw(k) * factor
       endif
!
! rain
!
       value = max(qmin,ncr(k,i,3))
       source = (-nraut(k)+nrcol(k)+niacr(k)+nsacr(k)+ngacr(k))*dtcld
       if (source>value) then
         factor = value/source
         nraut(k) = nraut(k) * factor
         nrcol(k) = nrcol(k) * factor
         niacr(k) = niacr(k) * factor
         nsacr(k) = nsacr(k) * factor
         ngacr(k) = ngacr(k) * factor
       endif
!
! update
!
       total =-(prevp(k)+psdep(k)+pgdep(k)+pigen(k)+pidep(k))
       q(k,i) = q(k,i)+total*dtcld
       qci(k,i,1) = max(qci(k,i,1) - (praut(k)+pracw(k) &
                   +paacw(k)+paacw(k))*dtcld,0.)
       qrs(k,i,1) = max(qrs(k,i,1) + (praut(k)+pracw(k) &
                   +prevp(k)-piacr(k)-pgacr(k)-psacr(k))*dtcld,0.)
       qci(k,i,2) = max(qci(k,i,2) - (psaut(k)+praci(k) &
                   +psaci(k)+pgaci(k)-pigen(k)-pidep(k))*dtcld,0.)
       qrs(k,i,2) = max(qrs(k,i,2) + (psdep(k)+psaut(k)+paacw(k) &
                   -pgaut(k)+piacr(k)*delta3 + praci(k)*delta3 &
                   +psaci(k)-pgacs(k)-pracs(k)*(1.-delta2) + psacr(k)*delta2) &
                   *dtcld,0.)
       qrs(k,i,3) = max(qrs(k,i,3) + (pgdep(k)+pgaut(k)+piacr(k)*(1.-delta3) &
                   +praci(k)*(1.-delta3) + psacr(k)*(1.-delta2) &
                   +pracs(k)*(1.-delta2) + pgaci(k)+paacw(k) &
                   +pgacr(k)+pgacs(k))*dtcld,0.)
       ncr(k,i,2) = max(ncr(k,i,2) + (-nraut(k)-nccol(k)-nracw(k) &
                   -naacw(k)-naacw(k))*dtcld,0.)
       ncr(k,i,3) = max(ncr(k,i,3) + (nraut(k)-nrcol(k)-niacr(k) &
                   -nsacr(k)-ngacr(k))*dtcld,0.)
       value = - xls*(psdep(k)+pgdep(k)+pidep(k)+pigen(k)) &
               - xlv(k,i)*prevp(k) - xlf(k,i)*(piacr(k)+paacw(k) &
               + paacw(k)+pgacr(k)+psacr(k))
       t(k,i) = t(k,i) - value/cpm(k,i)*dtcld
!
     else
!
! cloud water
!
       value = max(qmin,qci(k,i,1))
       source=(praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
       if (source>value) then
         factor = value/source
         praut(k) = praut(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
       endif
!
! rain
!
       value = max(qmin,qrs(k,i,1))
       source = (-paacw(k)-praut(k)+pseml(k)+pgeml(k)-pracw(k) &
                 -paacw(k)-prevp(k))*dtcld
       if (source>value) then
         factor = value/source
         praut(k) = praut(k) * factor
         prevp(k) = prevp(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
         pseml(k) = pseml(k) * factor
         pgeml(k) = pgeml(k) * factor
       endif
!
! snow
!
       value = max(qmin,qrs(k,i,2))
       source=(pgacs(k)-pseml(k)-psevp(k))*dtcld
       if (source>value) then
         factor = value/source
         pgacs(k) = pgacs(k) * factor
         psevp(k) = psevp(k) * factor
         pseml(k) = pseml(k) * factor
       endif
!
! graupel
!
       value = max(qmin,qrs(k,i,3))
       source=-(pgacs(k)+pgevp(k)+pgeml(k))*dtcld
       if (source>value) then
         factor = value/source
         pgacs(k) = pgacs(k) * factor
         pgevp(k) = pgevp(k) * factor
         pgeml(k) = pgeml(k) * factor
       endif
!
! cloud
!
       value = max(qmin,ncr(k,i,2))
       source = (+nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
       if (source>value) then
         factor = value/source
         nraut(k) = nraut(k) * factor
         nccol(k) = nccol(k) * factor
         nracw(k) = nracw(k) * factor
         naacw(k) = naacw(k) * factor
       endif
!
! rain
!
       value = max(qmin,ncr(k,i,3))
       source = (-nraut(k)+nrcol(k)-nseml(k)-ngeml(k))*dtcld
       if (source>value) then
         factor = value/source
         nraut(k) = nraut(k) * factor
         nrcol(k) = nrcol(k) * factor
         nseml(k) = nseml(k) * factor
         ngeml(k) = ngeml(k) * factor
       endif
!
! update
!
       total = -(prevp(k)+psevp(k)+pgevp(k))
       q(k,i) = q(k,i) + total*dtcld
       qci(k,i,1) = max(qci(k,i,1) - (praut(k)+pracw(k) &
                   +paacw(k)+paacw(k))*dtcld,0.)
       qrs(k,i,1) = max(qrs(k,i,1) + (praut(k)+pracw(k) &
                   +prevp(k)+paacw(k)+paacw(k)-pseml(k)-pgeml(k))*dtcld,0.)
       qrs(k,i,2) = max(qrs(k,i,2) + (psevp(k)-pgacs(k) &
                   +pseml(k))*dtcld,0.)
       qrs(k,i,3) = max(qrs(k,i,3) + (pgacs(k)+pgevp(k) &
                   +pgeml(k))*dtcld,0.)
       ncr(k,i,2) = max(ncr(k,i,2) + (-nraut(k)-nccol(k)-nracw(k) &
                   -naacw(k)-naacw(k))*dtcld,0.)
       ncr(k,i,3) = max(ncr(k,i,3) + (nraut(k)-nrcol(k)+nseml(k) &
                   +ngeml(k))*dtcld,0.)
       value = -xlv(k,i)*(prevp(k)+psevp(k)+pgevp(k)) &
               -xlf(k,i)*(pseml(k)+pgeml(k))
       t(k,i) = t(k,i) - value/cpm(k,i)*dtcld
     endif
   enddo
!===============================================================================
!
! ccn activaiton and condensation/evaporation of clouds
!
!===============================================================================
!
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   if(ktopqr>0.0) lqr = .true.
   if(lqr) then
   ktop = ktopqr
   call slope_rain(1,kdim,ktop,qrs(:,i,1),den(:,i),denfac(:,i),t(:,i), &
      ncr(:,i,3), &
      rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:))
   do k = ktop, kts, -1
     if(qrs(k,i,1)>0.) then
       diar(k) = rslope(k,1)*drcoeff
!----------------------------------------------------------------
! nrtoc: conversion from rain to cloud [lh a14]
! (nr->nc)
       if(diar(k)<=di82) then
         ncr(k,i,2) = ncr(k,i,2) + ncr(k,i,3)
         ncr(k,i,3) = 0.
!----------------------------------------------------------------
! prtoc: conversion from rain to cloud [lh a15]
! (qr->qc)
         qci(k,i,1) = qci(k,i,1) + qrs(k,i,1)
         qrs(k,i,1) = 0.
       endif
     endif
   enddo
   endif
!
   ktop = ktopini
   do k = ktop, kts, -1
     tr = ttp/t(k,i)
     logtr = log(tr)
     qsat_mul(k) = psat*exp(logtr*(xa)+xb*(1.-tr))
     qsat_mul(k) = min(qsat_mul(k),0.99*p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k) = max(q(k,i) / qsat_mul(k),qmin)
   enddo
!
   call find_cloud_top(1,kdim,ktopini,rh_mul(:), one_1,ktoprh)
!
   ktop = ktoprh
   do k = ktop, kts, -1
!---------------------------------------------------------------
! pcact: qv -> qc [lh 8] [kk 14]
! ncact: nccn -> nc [lh a2] [kk 12]
     if(rh_mul(k)>1.) then
       ncact(k) = max(0.,((ncr(k,i,1) + ncr(k,i,2)) &
                    *min(1.,((rh_mul(k)-1.)/satmax)**actk) - ncr(k,i,2)))*rdtcld
       ncact(k) =min(ncact(k),max(ncr(k,i,1),0.)*rdtcld)
       pcact(k) = min(4.*pi*denr*(actr*1.E-6)**3*ncact(k)/ &
                     (3.*den(k,i)),max(q(k,i),0.)*rdtcld)
       q(k,i) = max(q(k,i) - pcact(k)*dtcld,0.)
       qci(k,i,1) = max(qci(k,i,1) + pcact(k)*dtcld,0.)
       ncr(k,i,1) = max(ncr(k,i,1) - ncact(k)*dtcld,0.)
       ncr(k,i,2) = max(ncr(k,i,2) + ncact(k)*dtcld,0.)
       t(k,i) = t(k,i) + pcact(k)*xlv(k,i)/cpm(k,i)*dtcld
     endif
   enddo
!
!-------------------------------------------------------------------------------
   ktop = ktopini
   do k = ktop, kts, -1
     tr=ttp/t(k,i)
     qsat_mul(k)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
     qsat_mul(k) = min(qsat_mul(k),0.99*p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
   enddo
!-------------------------------------------------------------------------------
! pcond: condensational/evaporational rate of cloud water [hl a46]
!
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,rh_mul(:), one_1,ktoprh)
   ktop = max(ktopqc,ktoprh)
   do k = ktop, kts, -1
     value = ((max(q(k,i),qmin)-(qsat_mul(k))) /(1.+(xlv(k,i))*(xlv(k,i)) &
              /(rv*(cpm(k,i)))*(qsat_mul(k)) /((t(k,i))*(t(k,i)))))
     pcond(k) = min(max(value*rdtcld,0.),max(q(k,i),0.)*rdtcld)
     if(qci(k,i,1)>0. .and. value<0.) then
       pcond(k) = max(value,-qci(k,i,1))*rdtcld
     endif
!----------------------------------------------------------------
! ncevp: evpration of cloud number concentration [lh A3]
! (nc->nccn)
      if(pcond(k)==-qci(k,i,1)*rdtcld) then
        ncr(k,i,1) = ncr(k,i,1) + ncr(k,i,2)
        ncr(k,i,2) = 0.
      endif
     q(k,i) = q(k,i) - pcond(k)*dtcld
     qci(k,i,1) = max(qci(k,i,1)+pcond(k)*dtcld,0.)
     t(k,i) = t(k,i) + pcond(k)*xlv(k,i)/cpm(k,i)*dtcld
   enddo
!===============================================================================
!
! adjustment of nc and nr according to avaialble size of clouds and rain
!
!===============================================================================
   ktop = max(ktopqr,ktopqc)
   do k = ktop, kts, -1
     if(qrs(k,i,1)>=qrmin .and. ncr(k,i,3)>=nrmin) then
       value = exp(log(((pidnr*ncr(k,i,3)) &
                    /(den(k,i)*qrs(k,i,1))))*((.33333333)))
       if(value<=lamdarmin) then
         value = lamdarmin
         ncr(k,i,3) = den(k,i)*qrs(k,i,1)*value**3/pidnr
       elseif(value>=lamdarmax) then
         value = lamdarmax
         ncr(k,i,3) = den(k,i)*qrs(k,i,1)*value**3/pidnr
       endif
     endif
     if(qci(k,i,1)>=qmin .and. ncr(k,i,2)>=ncmin ) then
       value = exp(log(((pidnc*ncr(k,i,2)) &
                    /(den(k,i)*qci(k,i,1))))*((.33333333)))
       if(value<=lamdacmin) then
         value = lamdacmin
         ncr(k,i,2) = den(k,i)*qci(k,i,1)*value**3/pidnc
       elseif(value>=lamdacmax) then
         value = lamdacmax
         ncr(k,i,2) = den(k,i)*qci(k,i,1)*value**3/pidnc
       endif
     endif
   enddo
!-------------------------------------------------------------------------------
!
   enddo i_loop ! i-loops
!
!===============================================================================
!
   enddo inner_loop ! dtcldcr- loops
!
!===============================================================================
!
!-------------------------------------------------------------------------------
! compute reflectivity and effective radius for ufs
!
   if (mod(kdt,(int(fhswr)/int(delt)))==1) then
   do i = its, ite
     dbz(:) = 0.; re_qc(:) = 0.; re_qi(:) = 0.; re_qs(:) = 0.
     call ufs_mp_reflectivity (q(:,i), qrs(:,i,1), qrs(:,i,2), qrs(:,i,3), &
                         ncr(:,i,3), &
                         t(:,i), p(:,i), dbz, kts, kte, i, i)
     do k = kts, kte
       refl_10cm(i,k) = max(-35., dbz(k))
     enddo
!
     call ufs_mp_effective_radius (t(:,i), qci(:,i,1), qci(:,i,2), qrs(:,i,2), &
                         den(:,i), qmin, t0c, &
                         ncr(:,i,2), &
                         re_qc, re_qi, re_qs, kts, kte, i, i)
     do k = kts, kte
       re_cloud(i,k) = max(recmin, min(re_qc(k), recmax))*1.e6
       re_ice(i,k) = max(reimin, min(re_qi(k), reimax))*1.e6
       re_snow(i,k) = max(resmin, min(re_qs(k), resmax))*1.e6
     enddo
   enddo
   endif
!-----------------------------------------------------------------------------
! assign local to passing variables
!
   ktop = kte
   do k = kts, ktop
     do i = its, ite
       t1(i,k) = t(k,i)
       q1(i,k) = q(k,i)
       nn1(i,k) = ncr(k,i,1)
       nc1(i,k) = ncr(k,i,2)
       nr1(i,k) = ncr(k,i,3)
     enddo
   enddo
!
   do k = kts, ktop
     do i = its, ite
       qc1(i,k) = qci(k,i,1)
       qr1(i,k) = qrs(k,i,1)
       qi1(i,k) = qci(k,i,2)
       qs1(i,k) = qrs(k,i,2)
       qg1(i,k) = qrs(k,i,3)
     enddo
   enddo
!
!
   end subroutine mp_ufs_run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_rain(idim, kdim, ktop, qrs, den, denfac, t, &
                         ncr, &
                         rslope, rslopeb, rslope2, rslope3, vt)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: idim, kdim, ktop
   real(kind=kind_phys), dimension( idim , kdim), &
        intent(in ) :: &
                                                                          qrs, &
                                                                          ncr, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension( idim, kdim), &
        intent(inout ) :: &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
  real(kind=kind_phys) :: lamdar, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
   lamdar(x,y,z)= exp(log(((pidnr*z)/(x*y)))*((.33333333)))
!
   do i = 1,idim
     do k = 1,ktop
       if(qrs(i,k)<=qrmin .or. ncr(i,k)<=nrmin ) then
         rslope(i,k) = rslopermax
         rslopeb(i,k) = rsloperbmax
         rslope2(i,k) = rsloper2max
         rslope3(i,k) = rsloper3max
       else
         rslope(i,k) = min(1./lamdar(qrs(i,k),den(i,k),ncr(i,k)),1.e-3)
         rslopeb(i,k) = exp(log(rslope(i,k))*(bvtr))
         rslope2(i,k) = rslope(i,k)*rslope(i,k)
         rslope3(i,k) = rslope2(i,k)*rslope(i,k)
       endif
       vt(i,k) = pvtr*rslopeb(i,k)*denfac(i,k)
       if(qrs(i,k)<=0.0) vt(i,k) = 0.0
     enddo
   enddo
end subroutine slope_rain
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_graupel(idim, kdim, ktop, qrs, den, denfac, t, rslope, &
                         rslopeb, rslope2, rslope3, vt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer :: idim, kdim, ktop
   real(kind=kind_phys), dimension( idim , kdim), &
        intent(in ) :: &
                                                                          qrs, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension( idim, kdim), &
        intent(inout ) :: &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
  real(kind=kind_phys) :: lamdag, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
   lamdag(x,y)= sqrt(sqrt(pidn0g/(x*y))) ! (pidn0g/(x*y))**.25
!
   do i = 1,idim
     do k = 1,ktop
       if(qrs(i,k)<=qrmin)then
         rslope(i,k) = rslopegmax
         rslopeb(i,k) = rslopegbmax
         rslope2(i,k) = rslopeg2max
         rslope3(i,k) = rslopeg3max
       else
         rslope(i,k) = 1./lamdag(qrs(i,k),den(i,k))
         rslopeb(i,k) = exp(log(rslope(i,k))*(bvtg))
         rslope2(i,k) = rslope(i,k)*rslope(i,k)
         rslope3(i,k) = rslope2(i,k)*rslope(i,k)
       endif
       vt(i,k) = pvtg*rslopeb(i,k)*denfac(i,k)
       if(qrs(i,k)<=0.0) vt(i,k) = 0.0
     enddo
   enddo
end subroutine slope_graupel
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_cloud(idim, kdim, ktop, qrs, ncr, den, denfac, t, qmin, &
                          rslope, rslope2, rslope3)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer :: idim, kdim, ktop
   real(kind=kind_phys), dimension( idim , kdim), &
        intent(in ) :: &
                                                                          qrs, &
                                                                          ncr, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension( idim, kdim), &
        intent(inout ) :: &
                                                                       rslope, &
                                                                      rslope2, &
                                                                      rslope3
  real(kind=kind_phys) :: lamdac, x, y, z, supcol, qmin
  integer :: i, j, k
!----------------------------------------------------------------
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
   lamdac(x,y,z)= exp(log(((pidnc*z)/(x*y)))*((.33333333)))
!
   do i = 1,idim
     do k = 1,ktop
       if(qrs(i,k)<=qmin .or. ncr(i,k)<=ncmin )then
         rslope(i,k) = rslopecmax
         rslope2(i,k) = rslopec2max
         rslope3(i,k) = rslopec3max
       else
         rslope(i,k) = 1./lamdac(qrs(i,k),den(i,k),ncr(i,k))
         rslope2(i,k) = rslope(i,k)*rslope(i,k)
         rslope3(i,k) = rslope2(i,k)*rslope(i,k)
       endif
     enddo
   enddo
end subroutine slope_cloud
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine slope_snow(idim, kdim, ktop, qrs, den, denfac, t, rslope, rslopeb, &
                            rslope2, rslope3, vt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer :: idim, kdim, ktop
   real(kind=kind_phys), dimension( idim , kdim), &
         intent(in ) :: &
                                                                          qrs, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
   real(kind=kind_phys), dimension( idim , kdim), &
         intent(inout ) :: &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
   real(kind=kind_phys), parameter :: t0c = 273.15
   real(kind=kind_phys), dimension( idim , kdim ) :: &
                                                                       n0sfac
   real(kind=kind_phys) :: lamdas, x, y, z, supcol
   integer :: i, j, k
!----------------------------------------------------------------
! size distributions: (x=mixing ratio, y=air density):
! valid for mixing ratio > 1.e-9 kg/kg.
   lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y))) ! (pidn0s*z/(x*y))**.25
!
   do i = 1,idim
     do k = 1, ktop
       supcol = t0c-t(i,k)
!
! n0s: intercept parameter for snow [m-4] [hdc 6]
!
       n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
       if(qrs(i,k)<=qrmin)then
         rslope(i,k) = rslopesmax
         rslopeb(i,k) = rslopesbmax
         rslope2(i,k) = rslopes2max
         rslope3(i,k) = rslopes3max
       else
         rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
         rslopeb(i,k) = exp(log(rslope(i,k))*(bvts))
         rslope2(i,k) = rslope(i,k)*rslope(i,k)
         rslope3(i,k) = rslope2(i,k)*rslope(i,k)
       endif
       vt(i,k) = pvts*rslopeb(i,k)*denfac(i,k)
       if(qrs(i,k)<=0.0) vt(i,k) = 0.0
     enddo
   enddo
end subroutine slope_snow
!-------------------------------------------------------------------
!
!
!-------------------------------------------------------------------
   subroutine semi_lagrange_sedim(im, km, ktop1, dendl, denfacl, tkl, dzl, &
                                  wwl, rql, precip, dt, lat, lon)
!-------------------------------------------------------------------
!
! this routine is a semi-lagrangain forward advection for hydrometeors
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_lagrangian
!
! input :
! im, km - dimension
! dendl - dry air density
! denfacl - sqrt(den0/den)
! tkl - temperature in k
! dzl - depth of model layer in meter
! wwl - terminal velocity at model layer in m/s
! dt - time step
! ktop - top layer for computing
! id - kind of precip: 0 test case; 1 raindrop 2: snow
! iter - how many time to guess mean terminal velocity: 0 pure forward.
! 0 : use departure wind for advection
! 1 : use mean wind for advection
! > 1 : use mean wind after iter-1 iterations
!
! inout :
! precip - precipitation
! rql - cloud density*mixing ratio
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
! implemented by song-you hong
! reference: Juang, H.-M., and S.-Y. Hong, 2010: Forward semi-Lagrangian advection
! with mass conservation and positive definiteness for falling
! hydrometeors. Mon. Wea. Rev., 138, 1778-1791
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer , intent(in ) :: im, km, ktop1
   integer , intent(in ) :: lat, lon
   real(kind=kind_phys) , intent(in ) :: dt
   real(kind=kind_phys), dimension(im,km), intent(in ) :: dzl
   real(kind=kind_phys), dimension(im,km), intent(in ) :: dendl
   real(kind=kind_phys), dimension(im,km), intent(in ) :: denfacl
   real(kind=kind_phys), dimension(im,km), intent(in ) :: tkl
   real(kind=kind_phys), dimension(km), intent(in ) :: wwl
   real(kind=kind_phys), dimension(km), intent(inout) :: rql
   real(kind=kind_phys), intent(inout) :: precip
!
! local variables
!
   real(kind=kind_phys), dimension(km) :: &
                                                        dz, &
                                                        ww, &
                                                        qq, &
                                                        wd, &
                                                        wa, &
                                                       den, &
                                                    denfac, &
                                                        tk, &
                                                        qn
   real(kind=kind_phys), dimension(km+1) :: wi, &
                                                        zi, &
                                                       qmi, &
                                                       qpi, &
                                                       dza, &
                                                        qa
   real(kind=kind_phys), dimension(km+2) :: za
!
   integer :: ktop, i, k, n, m, kk, kb, kt, iter, lond, latd
   real(kind=kind_phys) :: &
                            tl, tl2, qql, dql, qqd ,&
                            th, th2, qqh, dqh ,&
                            zsum, qsum, dim, dip ,&
                            zsumt, qsumt, zsumb, qsumb ,&
                            allold, allnew, zz, dzamin, cflmax, decfl
   real(kind=kind_phys) :: tmp
!
   real(kind=kind_phys), parameter :: &
                            cfac = 0.05, fa1 = 9./16., fa2 = 1./16.
!
   lond = 101
   latd = 1
!
   zi(:) = 0. ; wi(:) = 0. ; qa(:) = 0.
   wa(:) = 0. ; wd(:) = 0. ; dza(:) = 0.
   precip = 0.0 ; tmp = 0.0
!
   ktop = max(ktop1,2)
!------------------------------------------------------------------
!
   semi_loop : do i = 1,im
!
! assign local variables
     dz(:) = dzl(i,:)
     den(:) = dendl(i,:)
     denfac(:) = denfacl(i,:)
     tk(:) = tkl(i,:)
     qq(:) = rql(:)
     ww(:) = wwl(:)
!
! skip for no precipitation for all layers
     allold = 0.0
     do k = 1,ktop
       allold = allold + qq(k)
     enddo
     if(allold<=0.0) then
       cycle semi_loop
     endif
!
! compute interface values
     zi(1)=0.0
     do k = 1,ktop
       zi(k+1) = zi(k)+dz(k)
     enddo
!
! save departure wind
     wd(:) = ww(:)
     n=1
     wi(1) = ww(1)
     wi(ktop+1) = ww(ktop)
     do k=2,ktop
       wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
     enddo
!
! 3rd order interpolation to get wi
     wi(1) = ww(1)
     wi(2) = 0.5*(ww(2)+ww(1))
     do k=3,ktop-1
       wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
     enddo
     wi(ktop) = 0.5*(ww(ktop)+ww(ktop-1))
     wi(ktop+1) = ww(ktop)
!
! terminate of top of raingroup
     do k=2,ktop
       if( ww(k)==0.0 ) wi(k)=ww(k-1)
     enddo
!
! diffusivity of wi
     do k=ktop,1,-1
       decfl = (wi(k+1)-wi(k))*dt/dz(k)
       if( decfl > cfac ) then
         wi(k) = wi(k+1) - cfac*dz(k)/dt
       endif
     enddo
!
! compute arrival point
     do k=1,ktop+1
       za(k) = zi(k) - wi(k)*dt
     enddo
     za(ktop+2) = zi(ktop+1) !hmhj
!
     do k=1,ktop+1 !hmhj
       dza(k) = za(k+1)-za(k)
     enddo
!
! compute deformation at arrival point
     do k=1,ktop
       tmp = qq(k)*dz(k)/dza(k)
       qa(k) = tmp
     enddo
     qa(ktop+1) = 0.0
!
! estimate values at arrival cell interface with monotone
     do k=2,ktop
       dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
       dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
       if( dip*dim<=0.0 ) then
         qmi(k)=qa(k)
         qpi(k)=qa(k)
       else
         qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
         qmi(k)=2.0*qa(k)-qpi(k)
         if( qpi(k)<0.0 .or. qmi(k)<0.0 ) then
           qpi(k) = qa(k)
           qmi(k) = qa(k)
         endif
       endif
     enddo
     qpi(1)=qa(1)
     qmi(1)=qa(1)
     qmi(ktop+1)=qa(ktop+1)
     qpi(ktop+1)=qa(ktop+1)
!
! interpolation to regular point
     qn = 0.0
     kb=1
     kt=1
     intp : do k=1,ktop
            kb=max(kb-1,1)
            kt=max(kt-1,1)
! find kb and kt
            if( zi(k)>=za(ktop+1) ) then
              exit intp
            else
              find_kb : do kk=kb,ktop
                        if( zi(k)<=za(kk+1) ) then
                          kb = kk
                          exit find_kb
                        else
                          cycle find_kb
                        endif
              enddo find_kb
              find_kt : do kk=kt,ktop+2 !hmhj
                        if( zi(k+1)<=za(kk) ) then
                          kt = kk
                          exit find_kt
                        else
                          cycle find_kt
                        endif
              enddo find_kt
              kt = kt - 1
!
! compute q with piecewise parabolic method
              if( kt==kb ) then
                tl=(zi(k)-za(kb))/dza(kb)
                th=(zi(k+1)-za(kb))/dza(kb)
                tl2=tl*tl
                th2=th*th
                qqd=0.5*(qpi(kb)-qmi(kb))
                qqh=qqd*th2+qmi(kb)*th
                qql=qqd*tl2+qmi(kb)*tl
               qn(k) = (qqh-qql)/(th-tl)
              else if( kt>kb ) then
                tl=(zi(k)-za(kb))/dza(kb)
                tl2=tl*tl
                qqd=0.5*(qpi(kb)-qmi(kb))
                qql=qqd*tl2+qmi(kb)*tl
                dql = qa(kb)-qql
                zsum = (1.-tl)*dza(kb)
                qsum = dql*dza(kb)
                if( kt-kb>1 ) then
                do m=kb+1,kt-1
                  zsum = zsum + dza(m)
                  qsum = qsum + qa(m) * dza(m)
                enddo
                endif
                th=(zi(k+1)-za(kt))/dza(kt)
                th2=th*th
                qqd=0.5*(qpi(kt)-qmi(kt))
                dqh=qqd*th2+qmi(kt)*th
                zsum = zsum + th*dza(kt)
                qsum = qsum + dqh*dza(kt)
                qn(k) = qsum/zsum
              endif
               cycle intp
             endif
!
     enddo intp
!
! rain out
     sum_precip: do k=1,ktop
                   if( za(k)<0.0 .and. za(k+1)<=0.0 ) then !hmhj
                      precip = precip + qa(k)*dza(k)
                      cycle sum_precip
                   else if ( za(k)<0.0 .and. za(k+1)>0.0 ) then !hmhj
                      th = (0.0-za(k))/dza(k) !hmhj
                      th2 = th*th !hmhj
                      qqd = 0.5*(qpi(k)-qmi(k)) !hmhj
                      qqh = qqd*th2+qmi(k)*th !hmhj
                      precip = precip + qqh*dza(k) !hmhj
                      exit sum_precip
                   endif
                   exit sum_precip
     enddo sum_precip
!
! replace the new values
     rql(:) = qn(:)
!
! ----------------------------------
   enddo semi_loop
!
end subroutine semi_lagrange_sedim
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine cldf_mps_diag(idim, kdim, t, p, q, qc, qi, dx,cldf, ktop)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in ) :: idim, kdim, ktop
   real(kind=kind_phys), dimension(idim), intent(in ) :: dx
   real(kind=kind_phys), dimension(idim, kdim), intent(in ) :: t, p, q
   real(kind=kind_phys), dimension(idim, kdim), intent(in ) :: qc
   real(kind=kind_phys), dimension(idim, kdim), intent(in ) :: qi
   real(kind=kind_phys), dimension(idim, kdim), intent(out ) :: cldf
   ! local
   integer :: i, k, kk
   real(kind=kind_phys), parameter :: cldmin = 50., cldmax = 100.
   real(kind=kind_phys), parameter :: clddiff = cldmax - cldmin
   real(kind=kind_phys), parameter :: cldf_min = 0.5
   real(kind=kind_phys) :: cv_w_min,cv_i_min,cvf_min
   real(kind=kind_phys) :: cv_w_max,cv_i_max,cvf_max
   real(kind=kind_phys) :: dxkm
!
   do k = 1,ktop
     do i = 1,idim
       cldf(i,k) = 0.
     enddo
   enddo
! diagnostic method (kiaps)
   do k = 1,ktop
     do i = 1, idim
       cv_w_min = 4.82*(max(0.,qc(i,k))*1000.)**0.94/1.04
       cv_w_max = 5.77*(max(0.,qc(i,k))*1000.)**1.07/1.04
       cv_i_min = 4.82*(max(0.,qi(i,k))*1000.)**0.94/0.96
       cv_i_max = 5.77*(max(0.,qi(i,k))*1000.)**1.07/0.96
       cvf_min = cv_i_min+cv_w_min
       cvf_max = cv_i_max+cv_w_max
       dxkm = dx(i)/1000.
       cldf(i,k)= ((dxkm-cldmin)*(cvf_max-cvf_min)+clddiff*cvf_min)/clddiff
     enddo
   enddo
   do k = 1,ktop
     do i = 1,idim
       cldf(i,k)=min(1.,max(0.,cldf(i,k)))
       if(qc(i,k)+qi(i,k)<1.e-6) then
         cldf(i,k) = 0.
       endif
       if(cldf(i,k)<0.01) cldf(i,k) = 0.
       if(cldf(i,k)>0.99) cldf(i,k) = 1.
       if(cldf(i,k)>=0.01 .and. cldf(i,k)<cldf_min) cldf(i,k) = cldf_min ! min
     enddo
   enddo
!
end subroutine cldf_mps_diag
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine ufsmpinit(den0,denr,dens,cl,cpv)
!-------------------------------------------------------------------------------
  implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys), intent(in) :: den0,denr,dens,cl,cpv
!
   pi = 4.*atan(1.)
   pisq = pi * pi
   xlv1 = cl-cpv
!
   pidnc = pi*denr/6.
!
   bvtr1 = 1.+bvtr
   bvtr2 = 2.+bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   bvtr6 = 6.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4) ! 17.837825
   g6pbr = rgmma(bvtr6)
   bvtr5 = 5.+bvtr
   bvtr7 = 7.+bvtr
   bvtr2o5 = 2.5+.5*bvtr
   bvtr3o5 = 3.5+.5*bvtr
   g2pbr = rgmma(bvtr2)
   g5pbr = rgmma(bvtr5)
   g7pbr = rgmma(bvtr7)
   g5pbro2 = rgmma(bvtr2o5)
   g7pbro2 = rgmma(bvtr3o5)
   pvtr = avtr*g5pbr/24.
   pvtrn = avtr*g2pbr
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   pidn0r = pi*denr*n0r
   pidnr = 4.*pi*denr
   precr1 = 2.*pi*1.56
   precr2 = 2.*pi*.31*avtr**.5*g7pbro2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1) !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4) ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r = pi*denr*n0r
   pidn0s = pi*dens*n0s
   pacrc = pi*n0s*avts*g3pbs*.25*eacrc
!
   bvtg1 = 1.+bvtg
   bvtg2 = 2.5+.5*bvtg
   bvtg3 = 3.+bvtg
   bvtg4 = 4.+bvtg
   g1pbg = rgmma(bvtg1)
   g3pbg = rgmma(bvtg3)
   g4pbg = rgmma(bvtg4)
   pacrg = pi*n0g*avtg*g3pbg*.25
   g5pbgo2 = rgmma(bvtg2)
   pvtg = avtg*g4pbg/6.
   precg1 = 2.*pi*n0g*.78
   precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
   pidn0g = pi*deng*n0g
!
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rslopegmax = 1./lamdagmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rslopegbmax = rslopegmax ** bvtg
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rslopeg2max = rslopegmax * rslopegmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
   rslopeg3max = rslopeg2max * rslopegmax
!
   rslopecmax = 1./lamdacmax
   rslopec2max = rslopecmax * rslopecmax
   rslopec3max = rslopec2max * rslopecmax
end subroutine ufsmpinit
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
real function rgmma(x)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
! rgmma function: use infinite product form
   real(kind=kind_phys) :: euler
   parameter (euler=0.577215664901532)
   real(kind=kind_phys) :: x, y
   integer :: i
   if(x==1.)then
     rgmma=0.
       else
     rgmma=x*exp(euler*x)
     do i=1,10000
       y=float(i)
       rgmma=rgmma*(1.000+x/y)*exp(-x/y)
     enddo
     rgmma=1./rgmma
   endif
end function rgmma
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine find_cloud_top(im,km,ktop,qq,value,ktopout)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer , intent(in ) :: im, km, ktop
   real(kind=kind_phys), dimension(km) , intent(in ) :: qq
   real(kind=kind_phys) , intent(in ) :: value
   integer , intent(inout) :: ktopout
   integer :: i,k
!
! do i = 1,im
    ktopout = 0
    find_qrtop : do k = ktop,1, -1
      if(qq(k)>value) then
        ktopout = k
        exit find_qrtop
      else
        cycle find_qrtop
      endif
    enddo find_qrtop
! enddo
  return
end subroutine find_cloud_top
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine ufs_mp_effective_radius (t, qc, qi, qs, rho, qmin, t0c, &
                                nc, &
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)
!-------------------------------------------------------------------------------
! compute radiation effective radii of cloud water, ice, and snow for
! single-moment microphysics.
! these are entirely consistent with microphysics assumptions, not
! constant or otherwise ad hoc as is internal to most radiation
! schemes.
! coded and implemented by soo ya bae, kiaps, january 2015.
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!..sub arguments
   integer, intent(in) :: kts, kte, ii, jj
   real(kind=kind_phys), intent(in) :: qmin
   real(kind=kind_phys), intent(in) :: t0c
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: t
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: qc
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: nc
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: qi
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: qs
   real(kind=kind_phys), dimension( kts:kte ), intent(in):: rho
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qc
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qi
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qs
!..local variables
   integer:: i,k
   integer :: inu_c
   real(kind=kind_phys), dimension( kts:kte ):: ni
   real(kind=kind_phys), dimension( kts:kte ):: rqc
   real(kind=kind_phys), dimension( kts:kte ):: rnc
   real(kind=kind_phys), dimension( kts:kte ):: rqi
   real(kind=kind_phys), dimension( kts:kte ):: rni
   real(kind=kind_phys), dimension( kts:kte ):: rqs
   real(kind=kind_phys) :: temp
   real(kind=kind_phys) :: lamdac
   real(kind=kind_phys) :: supcol, n0sfac, lamdas
   real(kind=kind_phys) :: diai ! diameter of ice in m
   real(kind=kind_phys) :: bfactor, bfactor2, bfactor3
   double precision :: lamc
   logical :: has_qc, has_qi, has_qs
!..minimum microphys values
   real(kind=kind_phys), parameter :: r1 = 1.e-12
   real(kind=kind_phys), parameter :: r2 = 1.e-6
!..mass power law relations: mass = am*d**bm
   real(kind=kind_phys), parameter :: bm_r = 3.0
   real(kind=kind_phys), parameter :: obmr = 1.0/bm_r
   real(kind=kind_phys), parameter :: nc0 = 3.e8
   real(kind=kind_phys), parameter :: rqi0 = 50.e-3 ! 50 g m-3
!-----------------------------------------------------------------------
   has_qc = .false.
   has_qi = .false.
   has_qs = .false.
   do k = kts, kte
     ! for cloud
     rqc(k) = max(r1, qc(k)*rho(k))
     rnc(k) = max(R2, nc(k)*rho(k))
     if (rqc(k)>R1 .and. rnc(k)>R2) has_qc = .true.
     if (rqc(k)>r1) has_qc = .true.
     ! for ice
     rqi(k) = max(r1, qi(k)*rho(k))
     temp = (rho(k)*max(qi(k),qmin))
     temp = sqrt(sqrt(temp*temp*temp))
     ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
     rni(k)= max(r2, ni(k)*rho(k))
     if (rqi(k)>r1 .and. rni(k)>r2) has_qi = .true.
     ! for snow
     rqs(k) = max(r1, qs(k)*rho(k))
     if (rqs(k)>r1) has_qs = .true.
   enddo
   if (has_qc) then
     do k=kts,kte
       if (rqc(k)<=R1 .or. rnc(k)<=R2) CYCLE
       lamc = (pidnc*nc(k)/rqc(k))**obmr
       re_qc(k) = max(recmin,min(0.5*(1./lamc),recmax))
     enddo
   endif
  if (has_qi) then
     do k=kts,kte
       if (rqi(k)<=r1 .or. rni(k)<=r2) cycle
       temp = t0c - t(k)
       bfactor = -2.0 + 1.0e-3*temp*sqrt(temp)*log10(rqi(k)/rqi0)
       bfactor2 = bfactor*bfactor
       bfactor3 = bfactor2*bfactor
       temp = 377.4 + 203.3*bfactor+ 37.91*bfactor2 + 2.3696*bfactor3
       re_qi(k) = max(reimin,min(temp*1.e-6,reimax))
     enddo
   endif
   if (has_qs) then
     do k=kts,kte
       if (rqs(k)<=r1) cycle
       supcol = t0c-t(k)
       n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
       lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
       re_qs(k) = max(resmin,min(0.5*(1./lamdas), resmax))
     enddo
   endif
end subroutine ufs_mp_effective_radius
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine ufs_mp_reflectivity (qv1d, qr1d, qs1d, qg1d, &
                          nr1d, &
            t1d, p1d, dbz, kts, kte, ii, jj)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!..sub arguments
   integer, intent(in):: kts, kte, ii, jj
   real(kind=kind_phys), dimension(kts:kte), intent(in):: t1d, p1d
   real(kind=kind_phys), dimension(kts:kte), intent(in):: qv1d, qr1d, qs1d, qg1d
   real(kind=kind_phys), dimension(kts:kte), intent(in):: nr1d
   real(kind=kind_phys), dimension(kts:kte), intent(inout):: dbz
!..local variables
   real(kind=kind_phys), dimension(kts:kte):: temp, pres, qv, rho
   real(kind=kind_phys), dimension(kts:kte):: rr, nr, rs, rg
   real(kind=kind_phys):: temp_c
!
   double precision, dimension(kts:kte):: ilamr, ilams, ilamg
   double precision, dimension(kts:kte):: n0_r, n0_s, n0_g
   double precision:: lamr, lams, lamg
   logical, dimension(kts:kte):: l_qr, l_qs, l_qg
!
   real(kind=kind_phys), dimension(kts:kte):: ze_rain, ze_snow, ze_graupel
   double precision:: fmelt_s, fmelt_g
!
   integer:: i, k, k_0, kbot, n
   logical:: melti
!
   double precision:: cback, x, eta, f_d
   real(kind=kind_phys), parameter:: r=287.
!
!+---+
!
   do k = kts, kte
     dbz(k) = -35.0
   enddo
!+---+-----------------------------------------------------------------+
!..put column of data into local arrays.
!+---+-----------------------------------------------------------------+
   do k = kts, kte
     temp(k) = t1d(k)
     temp_c = min(-0.001, temp(k)-273.15)
     qv(k) = max(1.e-10, qv1d(k))
     pres(k) = p1d(k)
     rho(k) = 0.622*pres(k)/(r*temp(k)*(qv(k)+0.622))
!
     if (qr1d(k) > 1.e-9) then
       rr(k) = qr1d(k)*rho(k)
       nr(k) = nr1d(k)*rho(k)
       lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
       N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
       ilamr(k) = 1./lamr
       l_qr(k) = .true.
     else
       rr(k) = 1.e-12
       nr(k) = 1.e-12
       l_qr(k) = .false.
     endif
!
     if (qs1d(k) > 1.e-9) then
       rs(k) = qs1d(k)*rho(k)
       n0_s(k) = min(n0smax, n0s*exp(-alpha*temp_c))
       lams = (xam_s*xcsg(3)*n0_s(k)/rs(k))**(1./xcse(1))
       ilams(k) = 1./lams
       l_qs(k) = .true.
     else
       rs(k) = 1.e-12
       l_qs(k) = .false.
     endif
!
     if (qg1d(k) > 1.e-9) then
       rg(k) = qg1d(k)*rho(k)
       n0_g(k) = n0g
       lamg = (xam_g*xcgg(3)*n0_g(k)/rg(k))**(1./xcge(1))
       ilamg(k) = 1./lamg
       l_qg(k) = .true.
     else
       rg(k) = 1.e-12
       l_qg(k) = .false.
     endif
   enddo
!
!+---+-----------------------------------------------------------------+
!..locate k-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
   melti = .false.
   k_0 = kts
   do k = kte-1, kts, -1
      if ( (temp(k)>273.15) .and. l_qr(k) .and. l_qs(k+1) ) then
         k_0 = max(k+1, k_0)
         melti=.true.
         goto 195
      endif
   enddo
 195 continue
!+---+-----------------------------------------------------------------+
!..assume rayleigh approximation at 10 cm wavelength. rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. integrations of m(d)*m(d)*n(d)*dd.
!+---+-----------------------------------------------------------------+
   do k = kts, kte
      ze_rain(k) = 1.e-22
      ze_snow(k) = 1.e-22
      if (l_qr(k)) ze_rain(k) = n0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
      if (l_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/pi)*(6.0/pi) &
                              * (xam_s/900.0)*(xam_s/900.0) &
                              * n0_s(k)*xcsg(4)*ilams(k)**xcse(4)
   enddo
!+---+-----------------------------------------------------------------+
!..special case of melting ice (snow/graupel) particles. assume the
!.. ice is surrounded by the liquid water. fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. uses code from uli blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+
   if (melti .and. k_0>=kts+1) then
     do k = k_0-1, kts, -1
!..reflectivity contributed by melting snow
       if (l_qs(k) .and. l_qs(k_0) ) then
        fmelt_s = max(0.005d0, min(1.0d0-rs(k)/rs(k_0), 0.99d0))
        eta = 0.d0
        lams = 1./ilams(k)
        do n = 1, nrbins
           x = xam_s * xxds(n)**xbm_s
           call rayleigh_soak_wetgraupel (x,dble(xocms),dble(xobms), &
                 fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                 cback, mixingrulestring_s, matrixstring_s, &
                 inclusionstring_s, hoststring_s, &
                 hostmatrixstring_s, hostinclusionstring_s)
           f_d = n0_s(k)*xxds(n)**xmu_s * dexp(-lams*xxds(n))
           eta = eta + f_d * cback * simpson(n) * xdts(n)
        enddo
        ze_snow(k) = sngl(lamda4 / (pi5 * k_w) * eta)
       endif
     enddo
   endif
   do k = kte, kts, -1
     dbz(k) = 10.*log10((ze_rain(k)+ze_snow(k))*1.d18)
   enddo
end subroutine ufs_mp_reflectivity
end module mp_ufs
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
