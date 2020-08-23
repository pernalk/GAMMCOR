      subroutine service_cvb
      implicit double precision (a-h,o-z)
      parameter (nstrin=4)
      character*8 keyword,string
      include "common/bas"
      include "common/tapes"
      include "common/maxatm"
      include "common/dumpinfow"
      include "common/big"
      dimension rinfo(7),string(nstrin)
      data string/'WRITE','START','DEBUG','SPECIAL'/
      write(iout,6100) nsk
      write(iout,6200) 'basis functions',ntg
C Defaults
      irec=2100
      iunit=2
      iwrite1=iout
      iwrite2=0
      iwrite3=0
      idim1=0
      idim2=0
      idim3=0
      idim4=0
      idebug=0
      call backspace_inp
      call input(ncol)
      call getx(1,keyword,1)
      call fzero(rinfo,7)
      do 7042 icol=2,ncol
7042  call inpf(icol,rinfo(icol-1))
C RINFO: irec.iunit,iwrite1,iwrite2,iwrite3 (see also 'start' keyword)
7002  info=int(rinfo(1))
      if(info.gt.0) irec=info
      info=int(10d0*(rinfo(1)-dble(info))+0.1d0)
      if(info.gt.0) iunit=info
      info=int(rinfo(2))
      if(info.ne.0) iwrite1=info
      info=int(rinfo(3))
      if(info.ne.0) iwrite2=info
      info=int(rinfo(4))
      if(info.ne.0) iwrite3=info
C Process keywords
7000  call input(ncol)
      call getx(1,keyword,1)
      call fzero(rinfo,7)
      do 7043 icol=2,ncol
7043  call inpf(icol,rinfo(icol-1))
      do 7050 i=1,nstrin
      if(keyword(1:5).eq.string(i)(1:5)) goto 7100
7050  continue
      call backspace_inp
      goto 7777
7100  goto(7001,7002,7003,7004),i
C...write
7001  info=int(rinfo(1))
      iwrite1=info
      info=int(rinfo(2))
      iwrite2=info
      info=int(rinfo(3))
      iwrite3=info
      goto 7000
C...debug
7003  info=int(rinfo(1))
      if(info.ne.0) idebug=info
      goto 7000
C...special
7004  idim1=int(rinfo(1))
      idim2=int(rinfo(2))
      idim3=int(rinfo(3))
      idim4=int(rinfo(4))
      len=idim1*idim2+idim3
      if(len.gt.0) then
        write(iout,6250) irec,iunit
        ibuf1=icorr(len+idim4)
        call lesw(q(ibuf1),len+idim4,iunit,irec,0)
        call dlcvout_cvb(q(ibuf1+idim4),len,iwrite1,'special')
        call corlsr(ibuf1)
        return
      endif
      goto 7000
C
7777  if(iwrite1.eq.0 .and. iwrite2.eq.0 .and. iwrite3.eq.0)
     >  call caserr('nothing to write in service!')
      if(idebug.gt.0) write(iout,6250) irec,iunit
      if(nsk.ne.1) write(iout,6210) (nt(ii),ii=1,nsk)
      ibuf1=icorrz(ntg*ntg)
      ibuf2=icorr(ntqg)
      call read_info(irec,iunit,0,idum,keyword)
      iset=1
      do k=1,dump_nset
        if (dump_name(k).eq.'ORBITALS') iset=dump_set(k)
      end do
      write(6,*) 'Occupation numbers of',ntg,'orbitals'
      call read_occ(q(ibuf2),iset)
      write(6,'(6f12.8)') (q(ibuf2-1+ii),ii=1,ntg)
      write(iout,6250) irec,iunit,iset
      call read_orb(q(ibuf2),iset)
      do 2100 isk=1,nsk
      iget=ibuf2 + ntq(isk)
      do 2100 jbas=1,nt(isk)
      iput=ibuf1 + (jbas+nts(isk)-1)*ntg + nts(isk)
      call fmove(q(iget),q(iput),nt(isk))
2100  iget=iget+nt(isk)
      call corlsr(ibuf2)
      call dlcvout_cvb(q(ibuf1),ntg*ntg,iwrite1,
     >             'symmetry orbital basis')
      if(iwrite2.eq.0 .and. iwrite3.eq.0) then
        call corlsr(ibuf1)
        return
      endif
      call dlcaoso_cvb(nbfao,nbfso,ibuft,iwrite3,idebug)
      if(iwrite2.eq.0) then
        call corlsr(ibuf1)
        return
      endif
      if(nbfao.ne.nbfso) write(iout,6200) 'atomic orbitals',nbfao
      if(nbfao.eq.nbfso) call dlcchk_cvb(q(ibuft),nbfao,nbfso)
      ibuf2=icorrz(nbfao*nbfao)
      do 3100 imo=1,ntg
      itrn=ibuft-1
      iput=ibuf2-1 + (imo-1)*nbfao
      iget=ibuf1-1 + (imo-1)*nbfso
      do 3100 ibfso=1,nbfso
      do 3100 ibfao=1,nbfao
      itrn=itrn+1
3100  q(iput+ibfao)=q(iput+ibfao) + q(itrn)*q(iget+ibfso)
      call dlcvout_cvb(q(ibuf2),nbfao*nbfao,iwrite2,
     >             'atomic orbital basis')
      call corlsr(ibuf1)
      return
6100  format(' Number of symmetry types: ',i1)
6200  format(' Total number of ',a,':',i4)
6210  format(' Number of basis functions (by symmetry):',/,8i4)
6250  format(/,' Record',i5,' on MOLPRO file ',i1,:,
     >         '    (Orbital set',i2,')'/)
      end
      subroutine dlcvout_cvb(v,n,iwritev,text)
      implicit double precision (a-h,o-z)
      include "common/tapes"
      character*(*) text
      dimension v(n)
      if(iwritev.gt.0) then
        write(iwritev,'(4d20.12)') v
        write(iout,6100) text,' ',iwritev
      elseif(iwritev.lt.0) then
        write(-iwritev) v
        write(iout,6100) text,' as a binary file ',-iwritev
      endif
      return
      entry dlcvout2_cvb(iwritev,n,v)
      if(iwritev.gt.0) then
        write(iwritev,*) n
        write(iwritev,'(4d20.12)') v
      elseif(iwritev.lt.0) then
        write(-iwritev) n
        write(-iwritev) v
      endif
      return
6100  format(' Vectors (',a,') written',a,'to channel ',i2)
      end
      subroutine dlcchk_cvb(t,na,ns)
      implicit double precision (a-h,o-z)
      include "common/tapes"
      dimension t(na,ns)
      big=0.d0
      do 1100 j=1,na
      big=dmax1(dabs(ddot_X(ns,t(j,1),na,t(j,1),na)-1.d0),big)
      do 1100 i=1,j-1
1100  big=dmax1(dabs(ddot_X(ns,t(i,1),na,t(j,1),na)),big)
      if(big.gt.1.d-14) write(iout,*)
     >  'WARNING: Transposition reliability is',big
      return
      end
      subroutine dlcaoso_cvb(nbfao,nbfso,ibuft,iwrite3,idebug)
      implicit double precision (a-h,o-z)
      include "common/big"
      include "common/tapes"
      name=0
      ifil=1
      call train(nnmx,nbfao,nbfso,ip,ipp,ifil,name)
      if(nbfao.ne.nbfso) write(iout,6200) 'atomic orbitals',nbfao
      if(idebug.gt.0) then
        write(iout,*) 'AO-SO TRANSFORMATION:'
        write(iout,6400) (iq(ip+nn),q(ipp+nn),nn=0,nnmx-1)
      endif
      ibuft=icorrz(nbfao*nbfso)
      do 91100 nn=0,nnmx-1
91100 q(ibuft-1+iq(ip+nn))=q(ipp+nn)
      call dlcvout_cvb(q(ibuft),nbfso*nbfao,iwrite3,
     >             'AO-SO transformation')
      return
6200  format(' Total number of ',a,':',i4,/)
6400  format(4(i8,f12.5))
      end
      subroutine moscow_cvb
      implicit double precision (a-h,o-z)
      integer eta,contrc
      character(len=4) grpp
      include "common/big"
c
c.....read basis info
c
      nbasis=0
      call train(nnmx,nbfao,nbfso,itradd,itrmatd,1,nbasis)
      call basin(grpp,natoms,nuc,nprim,eta,itypf,iatmf,
     1     ncont,ncontc,nstart,nstop,ncct,mconu,maxseg,contrc,3,
     2     1,nbasis)
c
c.....modify addressing & normalization in transformation matrix
c
      call dfttra(iq(itypf),
     1  iq(itradd),q(itrmatd),iq(nstart),iq(ncct),ncont,nbfao,nnmx,
     2  nprim,ip2sr,ip2si,np2s,ip2sp,np2sp,iq(nstop),q(contrc),nbfso,
     3  q(eta),.true.)
c
c.....try to write MOSCOW input
c
      lquant=icori(nprim)
      icen=icori(nprim)
      call moscow2_cvb(iq(lquant),iq(icen),iq(itypf),
     > q(eta),q(nuc),q(ip2sr),iq(ip2si),nbfso,natoms,nprim,np2s)
      call corlsi(itradd)
      return
      end
      subroutine moscow2_cvb(lquant,icen,itypf,
     > eta,vlist,p2s,ip2s,nbfso,natoms,nprim,np2s)
      implicit double precision (a-h,o-z)
      parameter (lmax=4)
      include "common/big"
      include "common/tapes"
      dimension lquant(nbfso),icen(nbfso),itypf(nprim),
     > eta(nprim,4),vlist(natoms,4),p2s(np2s),ip2s(np2s)
      dimension lfind(0:lmax+1)
      data lfind/1,2,5,11,21,35/,cut/1.d-10/
      ierr=0
      ncontr=0
      do 1100 i=1,nprim
      if(itypf(i).gt.lfind(lmax+1)) then
        write(iout,*) 'Error: LQUANT too large for MOSCOW'
        ierr=ierr+1
      else
        do 1200 lq=0,lmax
        if(itypf(i).eq.lfind(lq)) then
          ncontr=ncontr+1
          lquant(ncontr)=lq
          nfunlq=(lq+1)*(lq+2)/2
          do 1300 j=1,nfunlq-1
          if(itypf(i+j).ne.itypf(i)+j) then
            write(iout,*) 'Error: ITYPE: ',itypf(i+j),itypf(i)+j
            ierr=ierr+1
          endif
          do 1400 k=1,4
          if(dabs(eta(i+j,k)-eta(i,k)).gt.cut) then
            write(iout,*) 'Error: ordering: ',eta(i,k),eta(i+j,k)
            ierr=ierr+1
          endif
1400      continue
1300      continue
          icen(ncontr)=0
          do 1500 j=1,natoms
          do 1600 k=1,3
1600      if(dabs(vlist(j,k)-eta(i,k)).gt.cut) goto 1500
          icen(ncontr)=j
1500      continue
          if(icen(ncontr).eq.0) then
            write(iout,*) 'Error: location: ',(eta(i,k),k=1,3)
            ierr=ierr+1
          endif
        endif
1200    continue
      endif
1100  continue

      if(ierr.eq.0) then
        write(iout,6100) ncontr,natoms,nbfso,nbfso
        write(iout,6200) (lquant(i),i=1,ncontr)
        write(iout,6300) (icen(i),i=1,ncontr)
        write(iout,6400) ncontr
        write(iout,6500) ((vlist(i,j),j=1,4),i=1,natoms)
        do 2100 i=1,nprim
        do 2100 lq=0,lmax
2100    if(itypf(i).eq.lfind(lq)) write(iout,6600) 1.d0,eta(i,4)
        write(iout,*) np2s
        do 2200 i=1,np2s
2200    write(iout,*) p2s(i),ip2s(i)
c        ibuft=icorrz(nprim*nbfso)
c        do 2300 i=1,np2s
c2300    q(ibuft+ip2s(i))=p2s(i)
c        call outsqr(q(ibuft),nprim,nprim,nbfso,'transformation')
         write(iout,*) 'All done!'
      else
        write(iout,*) 'Abandonned due to errors:',ierr
      endif
      return
6100  format(' &MOSCINP  NCONTR=',I3,', NATOMS=',I3,', MBAS=',I3,
     > ', ISP=+1,',/,' NELEC=0, MACT=',I3,', ITYPE= 1,')
6200  format(' LQUANT=',25(I1,','),255(:,/,8X,30(I1,',')))
6300  format(' ICEN=',15(I3,','),255(:,/,6X,15(I3,',')))
6400  format(' ICONTR=',I3,'*1,',/,' &END')
6500  format(1x,3f15.9,2x,f6.2)
6600  format(2f14.6)
      end
      subroutine skmdlc_cvb
      implicit double precision (a-h,o-z)
      parameter (nrex=8,nstrin=3,lbiel=1000)
      character*8 keyword,string,origin
      character*2 chrex(nrex)
      logical eij,ekl,phase
      dimension rinfo(7),string(nstrin),mrex(nrex)
      dimension ekern(1)
      include "common/bas"
      include "common/big"
      include "common/corbdim"
      include "common/maxbfn"
      include "common/molen"
      include "common/syminf"
      include "common/tapes"
      data zero/0.0d0/,half/0.5d0/,one/1.0d0/
      data mrex/1100,1200,6*1700/,phase/.false./
      data chrex/' S',' F',' X',' Y',' Z','LX','LY','LZ'/
      data string/'DEBUG','START','WRITE'/

      write(iout,6100) nsk
      write(iout,6200) 'basis functions',ntg
      if(nsk.ne.1) write(iout,6210) (nt(ii),ii=1,nsk)
      nint=ntg*(ntg+1)/2
      write(iout,6800) nint*(nint+1)/2,nint

c Defaults
      iunit1=1
      iunit2=1
      idebug=0
      ioref1=44
      ioref2=44
      iskip=0
      origin=' '
C Process keywords
7000  call input(ncol)
      call getx(1,keyword,1)
      call fzero(rinfo,7)
      do 7043 icol=2,ncol
7043  call inpf(icol,rinfo(icol-1))
      do 7050 i=1,nstrin
      if(keyword(1:5).eq.string(i)(1:5)) goto 7100
7050  continue
      call backspace_inp
      goto 7777
7100  goto(7001,7002,7003),i
C...debug
7001  info=int(rinfo(1))
      if(info.ne.0) idebug=info
      goto 7000
C...start
7002  info=int(rinfo(1))
      iunit1=info
      info=int(rinfo(2))
      iunit2=info
      goto 7000
C...write
7003  info=int(rinfo(1))
      ioref1=info
      info=int(rinfo(2))
      ioref2=info
      goto 7000
C
7777  if(iunit1.le.0 .or. ioref1.eq.0) goto 9200
      write(iout,6600) 'one-electron',iunit1
      ibuf3=icorr(nint)
      ibuf4=icorr(ntdg)
      call dlcvout2_cvb(ioref1,1,ekern)
      if(idebug.ne.0) write(iout,6901) ekern(1)
      do 9300 irex=1,nrex
      if(irex.eq.3) then
        phase=.true.
        if(idebug.ne.0) write(iout,6902) 'dcore'
        call dlcvout2_cvb(ioref1,3,dipk)
        if(idebug.ne.0) write(iout,6904) dipk
      endif
      if(idebug.ne.0) write(iout,6902) chrex(irex)
      if(irex.lt.3) then
        call les(q(ibuf4),ntdg,iunit1,mrex(irex))
        isyop=1
      elseif(iskip.lt.2) then
        if(irex.gt.5) then
          phase=.false.
          call read_op(q(ibuf4),opnuc,chrex(irex)//origin,
     +                 0,0,1,dum,isyop)
        else
          call read_op(q(ibuf4),opnuc,'DM',irex-2,0,1,dum,isyop)
        endif
      else
        isyop=0
        phase=.false.
      endif
      call fzero(q(ibuf3),nint)
      if(isyop.eq.1) then
        iput=ibuf3
        do 5100 isk=1,nsk
        iget=ibuf4+ntd(isk)
        do 5100 i=1,nt(isk)
        iput=iput+nts(isk)
        call fmove(q(iget),q(iput),i)
        iget=iget+i
5100    iput=iput+i
      elseif(isyop.gt.1) then
        iget=ibuf4
        do 5200 isk=1,nsk
        do 5200 jsk=isk+1,nsk
        if(mult(isk,jsk).eq.isyop) then
          do 5300 i=1,nt(isk)
          nti=ibuf3-1 + nts(isk)+i
          do 5300 j=1,nt(jsk)
          ntj=nts(jsk)+j
          q(ntj*(ntj-1)/2 + nti)=q(iget)
5300      iget=iget+1
        endif
5200    continue
      endif
      if(phase) call dscal_X(nint,-one,q(ibuf3),1)
      call dlcvout2_cvb(ioref1,nint,q(ibuf3))
9300  if(idebug.ne.0) write(iout,6904) (q(ibuf3+jj),jj=0,nint-1)
      call corlsr(ibuf4)
      call corlsr(ibuf3)
      write(iout,6650)

9200  if(iskip.gt.0 .or. iunit2.le.0 .or. ioref2.eq.0) return
      mput=1
      mint=0

      nbuf=ntr
      call AO_Integral_Matrix_Get_init(nbuf)
      ibuf1=icorr(nbuf)
      ibuf2=icorr(lbiel)
      ibuf3=icori(lbiel*4)
      write(iout,6600) 'sorted two-electron',iunit2

      do 2100 isk=1,nsk
      do 2200 i=1,nt(isk)
      do 2200 j=1,i
      eij=(i.eq.j)
      if(idebug.gt.1) write(iout,*) 'Call AO_Integral_Matrix_Get_All'
      call AO_Integral_Matrix_Get_All(i,j,isk,isk,isk,q(ibuf1),nbuf,
     +                                iaddr,jproc)
      if(idebug.gt.1) write(iout,*) 'Returned-1:',iaddr,jproc
      if(iaddr.ne.0) then
        ii=ibuf1+iaddr-1
        do 2300 k=1,i
        le=k
        if(k.eq.i) le=j
        do 2300 l=1,le
        value=q(ii)
        if(idebug.gt.2) then
          iposq=ii+1-ibuf1
          write(iout,6500) i,j,k,l,value,iposq,iposq-iaddr+1
        endif
        if(value.ne.zero) then
          mij=i*(i-1)/2+j
          mkl=k*(k-1)/2+l
          if(eij) value=value+value
          if(k.eq.l) value=value+value
          if(mij.eq.mkl) value=value+value
          ix=i+nts(isk)
          jx=j+nts(isk)
          kx=k+nts(isk)
          lx=l+nts(isk)
          call savintx_cvb(ix,jx,kx,lx,value,q(ibuf2),iq(ibuf3),lbiel,
     +                     mput,mint,idebug,ioref2)
        endif
2300    ii=ii+1
      endif
2200  continue
2100  continue

      do 3100 isk=2,nsk
      do 3100 jsk=1,isk-1
      do 3200 i=1,nt(isk)
      do 3200 j=1,i
      eij=(i.eq.j)
      if(idebug.gt.1) write(iout,*) 'Call AO_Integral_Matrix_Get_All'
      call AO_Integral_Matrix_Get_All(i,j,isk,isk,jsk,q(ibuf1),nbuf,
     +                                iaddr,jproc)
      if(idebug.gt.1) write(iout,*) 'Returned-2:',iaddr,jproc
      if(iaddr.ne.0) then
        ii=ibuf1+iaddr-1
        do 3300 k=1,nt(jsk)
        do 3300 l=1,k
        ekl=(k.eq.l)
        value=q(ii)
        aaaa=q(ii+1)
        bbbb=q(ii+2)
        if(idebug.gt.2) then
          iposq=ii+1-ibuf1
          write(iout,6500) i,j,k,l,value,iposq,iposq-iaddr+1
          write(iout,6500) i,k,j,l,aaaa,iposq+1,iposq-iaddr+2
          write(iout,6500) i,l,j,k,bbbb,iposq+2,iposq-iaddr+3
        endif
        ix=i+nts(isk)
        jx=j+nts(isk)
        kx=k+nts(jsk)
        lx=l+nts(jsk)
        if(eij) then
          value=value+value
          aaaa=aaaa+aaaa
        endif
        if(ekl) then
          value=value+value
          aaaa=aaaa+aaaa
        endif
        call savintx_cvb(ix,jx,kx,lx,value,q(ibuf2),
     +                   iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
        call savintx_cvb(ix,kx,jx,lx,(aaaa+bbbb)*half,q(ibuf2),
     +                   iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
        if(.not.(eij.or.ekl))
     +    call savintx_cvb(ix,lx,jx,kx,(aaaa-bbbb)*half,q(ibuf2),
     +                     iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
3300    ii=ii+3
      endif
3200  continue
3100  continue

      do 4100 iska=4,nsk
      do 4100 iskb=3,iska-1
      iskab=mult(iska,iskb)
      do 4200 iskc=2,iskb-1
      iskd=mult(iskab,iskc)
      if(iskd.ge.iskc) goto 4200

      do 4250 i=1,nt(iska)
      do 4250 j=1,nt(iskb)
      if(idebug.gt.1) write(iout,*) 'Call AO_Integral_Matrix_Get_All'
      call AO_Integral_Matrix_Get_All(i,j,iska,iskb,iskc,q(ibuf1),nbuf,
     +                                iaddr,jproc)
      if(idebug.gt.1) write(iout,*) 'Returned-3:',iaddr,jproc
      if(iaddr.ne.0) then
        ii=ibuf1+iaddr-1
        do 4300 k=1,nt(iskc)
        do 4300 l=1,nt(iskd)
        value=q(ii)
        aaaa=q(ii+1)
        bbbb=q(ii+2)
        if(idebug.gt.2) then
          iposq=ii+1-ibuf1
          write(iout,6500) i,j,k,l,value,iposq,iposq-iaddr+1
          write(iout,6500) i,k,j,l,aaaa,iposq+1,iposq-iaddr+2
          write(iout,6500) i,l,j,k,bbbb,iposq+2,iposq-iaddr+3
        endif
        ix=i+nts(iska)
        jx=j+nts(iskb)
        kx=k+nts(iskc)
        lx=l+nts(iskd)
        call savintx_cvb(ix,jx,kx,lx,value,q(ibuf2),
     +                  iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
        call savintx_cvb(ix,kx,jx,lx,(aaaa+bbbb)*half,q(ibuf2),
     +                  iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
        call savintx_cvb(ix,lx,jx,kx,(aaaa-bbbb)*half,q(ibuf2),
     +                  iq(ibuf3),lbiel,mput,mint,idebug,ioref2)
4300    ii=ii+3
      endif
4250  continue
4200  continue
4100  continue

      call savintx_cvb(0,0,0,0,half,q(ibuf2),iq(ibuf3),lbiel,
     +                 mput,mint,idebug,ioref2)

      call corlsi(ibuf3)
      call corlsr(ibuf2)
      call corlsr(ibuf1)
      write(iout,6300) mint-1,intrea()
      write(iout,6650)
      return
6100  format(' Number of symmetry types: ',i1)
6200  format(' Total number of ',a,':',i4,/)
6210  format(' Number of basis functions (by symmetry):',/,8i4,/)
6300  format(' Integrals=',i10,'; intrel= ',i1)
6500  format(1x,4i4,f15.8,' <---',i10,i8,i4)
6600  format(/' Processing ',a,' integrals from MOLPRO file',i2)
6650  format(' Processing completed')
6800  format(' jmax= ',i10,'; nint= ',i7)
6901  format(' ERN: ',g16.9)
6902  format(' ',a,' array:')
6904  format(1x,4g16.9)
      end
      subroutine savintx_cvb(ilab,jlab,klab,llab,value,buff,labels,
     +                       lbuff,mput,mint,idebug,ioref2)
      implicit double precision (a-h,o-z)
      dimension buff(lbuff),labels(4,lbuff)
      include "common/tapes"
      data zero/0.0d0/
      if(value.eq.zero) return
      mint=mint+1
      if(idebug.gt.1) write(iout,6500) ilab,jlab,klab,llab,value,mint
      labels(1,mput)=ilab
      labels(2,mput)=jlab
      labels(3,mput)=klab
      labels(4,mput)=llab
      buff(mput)=value
      if(mput.eq.lbuff) then
        if(ioref2.gt.0) then
          write(ioref2,*) lbuff
          write(ioref2,*) ((labels(i,j),i=1,4),buff(j),j=1,lbuff)
        elseif(ioref2.lt.0) then
          write(-ioref2) lbuff
          write(-ioref2) ((labels(i,j),i=1,4),buff(j),j=1,lbuff)
        endif
        mput=1
      elseif(ilab.eq.0) then
        if(ioref2.gt.0) then
          write(ioref2,*) mput-1
          write(ioref2,*) ((labels(i,j),i=1,4),buff(j),j=1,mput-1)
          write(ioref2,*) 0
        elseif(ioref2.lt.0) then
          write(-ioref2) mput-1
          write(-ioref2) ((labels(i,j),i=1,4),buff(j),j=1,mput-1)
          write(-ioref2) 0
        endif
        mput=1
      else
        mput=mput+1
      endif
      return
6500  format(1x,4i4,f15.8,' --->',i10)
      end
      subroutine dwrite_cvb(text,len,d)
      implicit double precision(a-h,o-z)
      character*(*) text
      dimension d(len)
      write(6,6101) text,len
      write(6,'(6d20.10)') d
      return
6101  format(/,1x,a,'-particle density matrix (',i8,' elements)',/)
      end
      subroutine dwrite_cvb_tofile(recnum,text,len,d,ioff)
      implicit double precision(a-h,o-z)
      character*(*) text
      dimension d(len)
      integer recnum,ioff
      call writem(d,len,2,recnum,ioff*len,text)
      return
6101  format(/,1x,a,'-particle density matrix (',i8,' elements)',/)
      end
      block data bdatad2gen_cvb
      include "common/d2gen_cvb"
      data d2gen/.true./
      end
