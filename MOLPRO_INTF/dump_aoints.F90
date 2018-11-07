module dumpintao

private 
public dump_aoints,dump_molpro_sapt

contains

subroutine dump_molpro_sapt(mon)
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      integer :: mon
      character(32) :: onefile,rdmfile,recname,orbname,str
      character(8) :: unit

      ! two-electron integrals 
      ibase=icorr(0)

      ! choose monomer
      if(mon==1) then
        onefile='AOONEINT_A'
        rdmfile='2RDMA'
        recname='CASORBA'
        orbname='MOLPRO_A.MOPUN'
      elseif(mon==2) then
        onefile='AOONEINT_B'
        rdmfile='2RDMB'
        recname='CASORBB'
        orbname='MOLPRO_B.MOPUN'
      endif     

      iS = icorr(ntdg)
      iT = icorr(ntdg)
      iV = icorr(ntdg)
      iH = icorr(ntdg)
      iOrbCas = icorr(ntqg)
      call readm(q(iS), ntdg, 1, 1100, 0, str)
      call readm(q(iT), ntdg, 1, 1400, 0, str)
      call readm(q(iV), ntdg, 1, 1410, 0, str)
      call readm(q(iH), ntdg, 1, 1200, 0, str)
      ! 1-el integrals
      call dump_oneints(trim(onefile),q(iS),q(iT),q(iV),q(iH))
      ! 1- and 2-RDM
      call dump_gamma(trim(rdmfile),2,7200,2,7100)
      ! orbitals
      call dump_orbs(q(iOrbCas),trim(recname),trim(orbname))

      if(mon==2) then 
        iex=iexcom_status()
        if(iex.eq.1) call excom(2)
        ! 2-el integrals
        call dump_twoints
        if (iex.eq.1) call excom(1)
      endif

      call corlsr(ibase)

end subroutine dump_molpro_sapt

subroutine dump_monomer(mon)
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      integer :: mon

      iS = icorr(ntdg)
      iT = icorr(ntdg)
      iV = icorr(ntdg)
      iH = icorr(ntdg)
      call readm(q(iS), ntdg, 1, 1100, 0, str)
      call readm(q(iT), ntdg, 1, 1400, 0, str)
      call readm(q(iV), ntdg, 1, 1410, 0, str)
      call readm(q(iH), ntdg, 1, 1200, 0, str)
      call dump_oneints('AOONEINT_A',q(iS),q(iT),q(iV),q(iH))

      call corlsr(iH)
      call corlsr(iV)
      call corlsr(iT)
      call corlsr(iS)

end subroutine dump_monomer

subroutine dump_aoints
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      character(32) :: str,namx
      character(8) :: unit

      ibase=icorr(0)
      iex=iexcom_status()
      if (iex.eq.1) call excom(2)

      call dump_twoints

      if (iex.eq.1) call excom(1)

      ! 1- and 2-RDM
      call dump_gamma('2RDM',2,7200,2,7100)
 
      ! one electron integrals
      iS = icorr(ntdg)
      iT = icorr(ntdg)
      iV = icorr(ntdg)
      iH = icorr(ntdg)
      call readm(q(iS), ntdg, 1, 1100, 0, str)
      call readm(q(iT), ntdg, 1, 1400, 0, str)
      call readm(q(iV), ntdg, 1, 1410, 0, str)
      call readm(q(iH), ntdg, 1, 1200, 0, str)
      ! add symmetries
      call dump_oneints('AOONEINT.mol',q(iS),q(iT),q(iV),q(iH))

      call corlsr(iH)
      call corlsr(iV)
      call corlsr(iT)
      call corlsr(iS)

!      ! orbitals
!
      iOrbCas = icorr(ntqg)
      call dump_orbs(q(iOrbCas),'CASORB','MOLPRO.MOPUN')

      call corlsr(ibase)

end subroutine dump_aoints

subroutine dump_basinfo(iunit,nstate)
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"

      integer :: iunit,nstate

      !print*, nocc,ncore,nclos,nval,nact
      !print*, 'Act:',iact(1:nsk)
      !print*, 'Closed:',iclos(1:nsk)
      !print*, 'Total:',nt(1:nsk)
      write(iunit) 'BASINFO '
      write(iunit) int(nsk,kind=4)
      write(iunit) int(nstate,kind=4)
      write(iunit) int(iclos(1:nsk),kind=4)
      write(iunit) int(iact(1:nsk),kind=4)
      write(iunit) int(nt(1:nsk),kind=4) 

end subroutine dump_basinfo 

subroutine dump_gamma(outfile,i2fil,i2recnum,i1fil,i1recnum)
      implicit double precision (a-h,o-z)
      include "common/maxatm"
      include "common/tapes"
!      include "common/dumpinfow"
!      include "common/cstate"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
!      include "common/big"
!      include "common/cbas"
!      include "common/clseg"
!      include "common/cref"
!      include "common/ctran2"
!      include "common/code"
!      include "common/cmpp"
!      include "common/d2gen_cvb"
      integer :: i2fil,i2recnum,i1fil,i1recnum
      character(*) :: outfile
      integer :: i1off,i2off,itr
      double precision :: d2mat(ic1d),d1mat(nact*nact)

      nact2 = nact*nact       

      call find_free_unit(ifil)
      open(unit=ifil,file=trim(outfile),form='unformatted')

      call dump_basinfo(ifil,nstate)

      write(ifil) '2RDM    '
      write(ifil) int(ic1d,kind=4),int(nstate,kind=4)
      i2off = 0
      do i=1,nstate
         itr = i*(i-1)/2 + i
         i2off = (itr-1)*ic1d
         call lesw(d2mat,ic1d,i2fil,i2recnum,i2off)
         !write(6,*) i,itr,i2off
         !write(6,*) d2mat(1:ic1d)
         write(ifil) int(i,kind=4) 
         write(ifil) d2mat(1:ic1d)
      enddo

      write(ifil) '1RDM    '
      write(ifil) int(nact,kind=4),int(nact2,kind=4),int(nstate,kind=4)
      i1off = 0
      do i=1,nstate
         itr = i*(i-1)/2 + i
         i1off = (itr-1)*nact2
         call lesw(d1mat,nact*nact,i1fil,i1recnum,i1off)
         !write(6,*) i,itr,i1off
         !write(6,*) d1mat(1:nact2)*0.5d0
         write(ifil) int(i,kind=4) 
         write(ifil) d1mat(1:nact2)
      enddo

      close(ifil)

end subroutine dump_gamma

subroutine dump_twoints
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      character(32) :: str,namx
      character(8) :: unit


      length=ntr
      call AO_Integral_Matrix_Get_init(length)
      ibuf=icorr(length)

      call find_free_unit(ifil)
      open(unit=ifil,file='AOTWOINT.mol',form='unformatted')

      write(ifil) int(nsk,kind=4)
      write(ifil) int(nt,kind=4)

      do iska=1,nsk
        if(nt(iska).eq.0) cycle 
        call dumpints1(ifil,q(ibuf),length,iska)
        !call printints1(q(ibuf),length,iska)
      enddo   

      intyp=nsk
      do iska=2,nsk
        do iskb=1,iska-1
          intyp=intyp+1
          if(nt(iska).eq.0) cycle 
          if(nt(iskb).eq.0) cycle 
          call dumpints2(ifil,q(ibuf),length,iska,iskb)
          !call printints2(q(ibuf),length,intyp,iska,iskb)
        enddo         
      enddo         

      do iska=4,nsk
        do iskb=3,iska-1
        iskab=mult(iska,iskb)
        do iskc=2,iskb-1
          iskd=mult(iskab,iskc)
          if(iskd.ge.iskc) cycle 
          intyp=intyp+1 
          if(nt(iska).eq.0) cycle  
          if(nt(iskb).eq.0) cycle
          if(nt(iskc).eq.0) cycle
          if(nt(iskd).eq.0) cycle
          call dumpints3(ifil,q(ibuf),length,iska,iskb,iskc,iskd)
          !call printints3(q(ibuf),length,intyp,iska,iskb,iskc,iskd)
        enddo
        enddo
      enddo

      close(ifil)

end subroutine dump_twoints

subroutine dump_oneints(infil,Smat,Tmat,Vmat,Hmat)
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      include "common/cgeom"
      character(*)     :: infil
      integer :: natm
      double precision :: enuc
      double precision :: Smat(*),Tmat(*),Vmat(*),Hmat(*)
      double precision,allocatable :: chrg(:),xyz(:,:) 

      enuc=calc_enuc()
      natm=calc_natom()
      allocate(chrg(natm),xyz(natm,3)) 
      ii=1
      do i=1,ncen
        if(charg(i).ne.0d0) then
           chrg(ii)=charg(i)
           xyz(ii,1:3)=rr(1:3,i)
           ii=ii+1
        endif
      enddo

      call find_free_unit(ifil)
      open(unit=ifil,file=infil,form='unformatted')
      write(ifil)
      write(ifil) int(nsk,kind=4),int(nt(1:nsk),kind=4),int(nts(1:nsk),kind=4)
      write(ifil) enuc
      write(ifil) 'OVERLAP '
      write(ifil) Smat(1:ntdg)
      write(ifil) 'KINETINT'
      write(ifil) Tmat(1:ntdg)
      write(ifil) 'POTENTAL'
      write(ifil) Vmat(1:ntdg)
      write(ifil) 'ONEHAMIL'
      write(ifil) Hmat(1:ntdg)
      write(ifil) 'ISORDK  '
      write(ifil) int(natm,kind=4)
      write(ifil) chrg(1:natm),xyz(1:natm,1:3)
     
      close(ifil)
     
      deallocate(chrg,xyz)

contains 
      function calc_enuc() result(enuc)

      double precision :: dist
      double precision :: enuc

      enuc=0d0
      do ib=1,ncen-1
         do ia=ib+1,ncen
            dist = sqrt(sum((rr(:,ia)-rr(:,ib))**2))
            enuc = enuc + charg(ia)*charg(ib)/dist
         enddo
      enddo

      end function calc_enuc
      
      function calc_natom() result(natm)
 
      integer :: i,natm 
     
      natm = 0 
      do i=1,ncen 
         if(charg(i)/=0d0) natm = natm + 1
      enddo

      end function calc_natom

end subroutine dump_oneints

subroutine dump_orbs(OrbCas,recname,outfile)
      implicit double precision (a-h,o-z)
      include "common/corb"
      include "common/cbas"
      include "common/tapes"
      include "common/big"
      double precision :: OrbCas(*)
      character(*)     :: recname,outfile
      double precision :: OccNat(ntg)
      double precision :: test(ntg,ntg)
      character(8) :: unit,method


      call find_rec(recname,icasrec,icasfil)

      call read_info(icasrec, icasfil, 0, idiffCAS, method)
      !write(6,*) idiffCAS,method
      call read_orb(OrbCas, 1)
      call flush_dump
       
!      call find_rec('NATORB',inatrec,inatfil)
!      call read_info(inatrec, inatfil, 0, idiffNAT, method)
!      call read_orb(OrbNat, 1)
!      call read_occ(OccNat,1)
!      print*,  OccNat(1:ntg)
!      call flush_dump

      call find_free_unit(ifil)
      open(unit=ifil,file=outfile,form='unformatted')
      write(ifil) 'CASORB  '
      write(ifil) int(nsk,kind=4),int(nt(1:nsk),kind=4),int(nts(1:nsk),kind=4)
      write(ifil) OrbCas(1:ntqg)
      !write(ifil) 'NATORB'
      !write(ifil) nsk,nt,nts
      !write(ifil) OrbNat(1:ntqd)
      !write(ifil) OccNat(1:ntg)
      close(ifil)

end subroutine dump_orbs 

subroutine find_rec(str,irec,ifil)
      implicit double precision (a-h,o-z)
      character(*) :: str
      integer :: irec,ifil
      character(8) :: unit
    
      ii1=1
      ityp=1  
      irec=0
      ifil=0
      xtmp=0
      call getvar(str,xtmp,unit,ityp,ii1,ii1,ii1)
      irec=int(xtmp)
      ifil=int((xtmp-dble(irec))*10.1d0)

      if(irec==0.and.trim(str)=='NATORB') then
        write(6,*) 'WARNING: NATORB label not defined! Assuming rec=2140.2!'
        irec=2140
        ifil=2
      elseif(irec==0.and.trim(str)=='CASORB') then
        write(6,*) 'WARNING: CASORB label not defined! Assuming rec=2144.2!'
        irec=2144
        ifil=2
      endif

end subroutine find_rec

!-----------------------------------------------------------------------
      subroutine dumpints1(iunit,t,length,iska)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(1)

      write(iunit) int(iska,kind=4)

      do i=1,nt(iska)
        do j=1,i
         call AO_Integral_Matrix_Get(i,j,iska,iska,iska,t,length,ii)
         if (ii.gt.0) then
           lcd=i*(i-1)/2+j
           jj=ii-1
           do k=1,i-1
             jj=jj+k
             t(jj)=2*t(jj) 
           enddo
           if(i==j) then
             t(ii+lcd-1)=2*t(ii+lcd-1)
             t(ii:ii+lcd-1)=2*t(ii:ii+lcd-1)
           endif
           t(ii+lcd-1)=2*t(ii+lcd-1)
           write(iunit) int(i,kind=4),int(j,kind=4)
           write(iunit) t(ii:ii+lcd-1)
         endif
        enddo
      enddo

      write(iunit) int(-1,kind=4), int(-1,kind=4)

      end subroutine dumpints1

!-----------------------------------------------------------------------
      subroutine printints1(t,length,iska)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(1)

      write(iout,10) iska
10    format(/' Integral block=',i1)
      ioffint=int(blockadr(iska))

      ijkl=0
      do i=1,nt(iska)
        do j=1,i
         call AO_Integral_Matrix_Get(i,j,iska,iska,iska,t,length,ii)
         if (ii.gt.0) then
          do k=1,i
            le=k
            if(i.eq.k) le=j
            do l=1,le
              ijkl=ijkl+1
              write(iout,20) ijkl,ioffint+ijkl,i,j,k,l,t(ii)
20            format(1x,2i8,3x,4i4,f15.6)
              ii=ii+1
             end do
            end do
           end if
        end do
      end do
      return
      end

!-----------------------------------------------------------------------
      subroutine dumpints2(iunit,t,length,iska,iskb)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(*)

      write(iunit) int(iska,kind=4),int(iskb,kind=4)

      lcd=(nt(iskb)*(nt(iskb)+1))/2

      do i=1,nt(iska)
        do j=1,i
          call AO_Integral_Matrix_Get(i,j,iska,iska,iskb,t,length,ii)
          if (ii.gt.0) then

            if(i==j) then
              do m=ii,ii+lcd*3,3
                t2=(t(m+1)+t(m+2))
                t3=(t(m+1)-t(m+2))
                t(m+1)=t2             !(ik|jl)
                t(m+2)=t3             !(il|jk)
              enddo
            else
              m=ii-3
              do k=1,nt(iskb)
                do l=1,k-1
                  m=m+3 
                  t2=(t(m+1)+t(m+2))
                  t3=(t(m+1)-t(m+2))
                  t(m+1)=t2*0.5d0     !(ik|jl)
                  t(m+2)=t3*0.5d0     !(il|jk)
                enddo
                m=m+3 
                t2=(t(m+1)+t(m+2))
                t3=(t(m+1)-t(m+2))
                t(m+1)=t2             !(ik|jl)
                t(m+2)=t3             !(il|jk)
              enddo
            endif

            jj=ii-3
            do k=1,nt(iskb) 
              jj=jj+3*k  
              t(jj)=2*t(jj)
            enddo
            if(i==j) then
              t(ii:ii+3*lcd-1:3)=2*t(ii:ii+3*lcd-1:3)
              jj=ii-3  
              do k=1,nt(iskb)
                jj=jj+3*k  
                t(jj+1)=2*t(jj+1)
                t(jj+2)=2*t(jj+2)
              enddo
            endif

            write(iunit) int(i,kind=4),int(j,kind=4)
            write(iunit) t(ii:ii+3*lcd-1)

          end if
        end do
      end do

      write(iunit) int(-1,kind=4), int(-1,kind=4)

      end subroutine dumpints2

!-----------------------------------------------------------------------
      subroutine printints2(t,length,intyp,iska,iskb)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(*)

      write(iout,10) intyp,iska,iskb
10    format(/' Integral block=',i2,'  isym=',i1,'  jsym=',i1)

      ioffint=int(blockadr(intyp))
      lcd=(nt(iskb)*(nt(iskb)+1))/2

      ijkl=1
      do i=1,nt(iska)
        do j=1,i
         call AO_Integral_Matrix_Get(i,j,iska,iska,iskb,t,length,ii)
         if (ii.gt.0) then
          do m=ii+1,ii+lcd*3,3
           t2=(t(m)+t(m+1))
           t3=(t(m)-t(m+1))
           t(m)=t2              !(ik|jl)
           t(m+1)=t3            !(ik|jk)
          end do
          do k=1,nt(iskb)
            do l=1,k
              write(iout,20) ijkl,ioffint+ijkl,i,j,k,l,t(ii),t(ii+1),t(ii+2)
20            format(1x,2i8,3x,4i4,3f15.6)
              ijkl=ijkl+3
              ii=ii+3
            end do
          end do
          end if
        end do
      end do
      return
      end

!-----------------------------------------------------------------------
      subroutine dumpints3(iunit,t,length,iska,iskb,iskc,iskd)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(*)

      write(iunit) int(iska,kind=4),int(iskb,kind=4),int(iskc,kind=4),int(iskd,kind=4)

      lcd=nt(iskc)*nt(iskd)

      do i=1,nt(iska)
        do j=1,nt(iskb)
          call AO_Integral_Matrix_Get(i,j,iska,iskb,iskc,t,length,ii)
          if (ii.gt.0) then

            do m=ii,ii+lcd*3,3
              t2=(t(m+1)+t(m+2))
              t3=(t(m+1)-t(m+2))
              t(m+1)=t2*0.5d0         !(ik|jl)
              t(m+2)=t3*0.5d0         !(il|jk)
            end do

            write(iunit) int(i,kind=4),int(j,kind=4)
            write(iunit) t(ii:ii+3*lcd-1)

          end if
        end do
      end do

      write(iunit) int(-1,kind=4), int(-1,kind=4)

      end subroutine dumpints3

!-----------------------------------------------------------------------
      subroutine printints3(t,length,intyp,iska,iskb,iskc,iskd)
!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/cbas"
      include "common/corb"
      include "common/tapes"
      include "common/ciaib"
      dimension t(*)

      write(iout,10) intyp,iska,iskb,iskc,iskd
10    format(/' Integral block=',i2,'  isym=',i1,'  jsym=',i1,'  ksym=',i1,'  lsym=',i1)

      ioffint=int(blockadr(intyp))

      ijkl=1
      do i=1,nt(iska)
        do j=1,nt(iskb)
         call AO_Integral_Matrix_Get(i,j,iska,iskb,iskc,t,length,ii)
         if (ii.gt.0) then
          do m=ii+1,ii+nt(iskc)*nt(iskd)*3,3
           t2=(t(m)+t(m+1))
           t3=(t(m)-t(m+1))
           t(m)=t2              !(ik|jl)
           t(m+1)=t3            !(il|jk)
          end do
          do k=1,nt(iskc)
            do l=1,nt(iskd)

              write(iout,20) ijkl,ioffint+ijkl,i,j,k,l,t(ii),t(ii+1),t(ii+2)
20            format(1x,2i8,3x,4i4,3f15.6)
              ijkl=ijkl+3
              ii=ii+3
            end do
          end do
          end if
        end do
      end do
      return
      end

end module dumpintao

