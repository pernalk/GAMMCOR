      subroutine mprop1
      USE outputResult
      implicit double precision(a-h,o-z)
      parameter (minpr=1,maxpr=6)
      character*20 chars(minpr:maxpr)
      character*16 type,str
      logical many,first,lexpec
      character*22 string(2)
      common/cjsym/ jstsym,isymd
      include "common/maxatm"
      include "common/tapes"
      include "common/dumpinfow"
      include "common/cstate"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      include "common/cbas"
      include "common/clseg"
      include "common/cref"
      include "common/ctran2"
      include "common/code"
      include "common/cmpp"
      include "common/d2gen_cvb"
      common/cener/ ener(mxstsy,mxstt)
      dimension dt(mxstt),pmat(mxstt,mxstt),dmat(mxstt,maxpr),
     1          iqpos(8),orig(3),erel(mxstt*mxstsy)
      data string/' ','(State averaged value)'/
      data chars /'Overlap'
     1           ,'Kinetic energy'
     2           ,'One electron energy'
     3           ,'X dipole'
     4           ,'Y dipole'
     5           ,'Z dipole'
     6/
      if(iguess.eq.0) call caserr('NO CI VECTOR PROVIDED')
      ibase = icorr(0)
      debye=get_scalar('TODEBYE')
      namo=0
      ioff=0
      do 2 isym=1,nirr
2     iqpos(isym) = icorr(nsymao(isym)**2)
      call qget (iqpos)
      write(iout,'()')
      idpos=icorr(ntdgcx(1))
      idpos_spin=icorr(0)
      call fzero(q(idpos),ntdgcx(1))
      if(dmsav.gt.1.1d0) then
        idpos_spin=icorr(ntdgcx(1))
        call fzero(q(idpos_spin),ntdgcx(1))
      end if
      igam= icorr(ne)
      icivec = icorr(lenciv)
      call clearvar('DMX')
      call clearvar('DMY')
      call clearvar('DMZ')
      call clearvar('ENERGY')
      many=nstsym.gt.1
      do 1 istsym=1,nstsym
1     if(nstats(istsym).gt.1) many=.true.
      iset=0
      nstate_av=0
      do 200 istsym=1,nstsym
      call blank_line(iout)
      if(istsym.gt.1) call blank_line(iout)
      nci=ncisym(istsym)
      nstate=nstats(istsym)
      call cget (q(icivec),istsym,nstate)
      ivec = icivec
      npri = nprim+nfreez
      nco = ncore+nfreez
      npr2 = 0
      do 10 i=1,nirr
10    npr2 = npr2 + (nactt(i)*(nactt(i)+1))/2
      iden = icorr(npr2 * (nstate*(nstate+1))/2 )
c
c...  form (transition) density matrices
      ii = iden
      do 50 istat1 = 1,nstate
      ivec1 = (istat1-1)*nci+ivec
      istate = istat1
      do 50 istat2 = 1,istat1
      ivec2 = (istat2-1)*nci+ivec
      if (dmsav.gt.1.1d0 .and. istat1.eq.istat2) then
c..    spin density
       q(igam)=0d0
       call denst1spin (q(ivec1),q(ivec2),q(igam))
cstart debug
c;       call outvec (q(igam+ic1d),ne-ic1d,'gam')
cend
       if (q(igam).ge.-100d0) call savdm(q(igam+ic1d),q(idpos_spin),
     &    q(iqpos(1)),iset,'SPIN',.true.)
      end if
      if(d2gen) then
        write(6,*) 'ioff=',ioff
        call denst2(q(ivec1),q(ivec2),q(igam))
        !call dwrite_cvb('2',ic1d,q(igam))
        call dwrite_cvb_tofile(7200,'2RDM',ic1d,q(igam),ioff)
        !call dwrite_cvb('1',nact*nact,q(igam+ic1d))
        call dwrite_cvb_tofile(7100,'1RDM',nact*nact,q(igam+ic1d),ioff)
        ioff = ioff + 1
c        d2gen=.false.
      endif
      call denst1 (q(ivec1),q(ivec2),q(igam))
cstart debug
c;       call outvec (q(igam+ic1d),ne-ic1d,'gam')
cend
      if(istat1.eq.istat2)
     1    call savdm(q(igam+ic1d),q(idpos),q(iqpos(1)),iset,'CHARGE',
     2    (.not.many.or.dmsav.gt.0d0))
c...  symmetry pack density matrix
      do 40 isym=1,nirr
      do 30 i=1,nact
      if (itypea(i).ne.isym) goto 30
      do 20 j=1,i
      if (itypea(j).ne.isym) goto 20
      q(ii) = q(igam+ic1d+(i-1)+(j-1)*nact)
      if (i.ne.j) q(ii) = q(ii) + q(igam+ic1d+(j-1)+(i-1)*nact)
      ii = ii + 1
20    continue
30    continue
40    continue
c      call outsqr (q(igam+ic1d),nact,nact,nact,'raw density')
c      call outvec (q(ii-npr2),npr2,'symmetrised density')
50    continue
c
c...  so now we can do properties
      iprop = icorr(npr2)
      call fzero (dt,nstate)
      call fzero(dmat,maxpr*mxstt)
      dtav = 0.0d0
      vt = -1.0d0
      do 160 i=1,n1elec
      ipr = i1elec(i)
      if (ipr.gt.maxpr .or. ipr.lt.minpr) then
      write(iout,70) ipr
70    format (/' *** REQUESTED PROPERTY NUMBER',i8,' IS OUT OF RANGE')
      goto 160
      end if
      call pget (q(iprop),ipr,znuc,zcor)
c     call outvec (q(iprop),npr2,'property integrals')
      ii = iden
      sign=1.0d0
      if(ipr.ge.4.and.ipr.le.6) sign=-1.0d0
      ityp=0
      nonz=0
      do 90 istate=1,nstate
      if(weight(istate,istsym).ne.0d0) ityp=ityp+1
      do 80 jstate=1,istate
      pmat(jstate,istate) = sign*ddot_X(npr2,q(ii),1,q(iprop),1)
      pmat(istate,jstate) = pmat(jstate,istate)
      if(istate.eq.jstate) goto 80
      if(pmat(istate,jstate).ne.0d0) nonz=nonz+1
80    ii = ii + npr2
      pmat(istate,istate) = pmat(istate,istate) + znuc + zcor*sign
      if(pmat(istate,istate).ne.0d0) nonz=nonz+1
90    dmat(istate,ipr) = pmat(istate,istate)
c...  dipole
      if (ipr.ge.4.and.ipr.le.6) then
      do 100 istate=1,nstate
100   dt(istate)=dt(istate)+dmat(istate,ipr)**2
      end if
c..   kinetic energy
      pav=ddot_X(nstate,dmat(1,ipr),1,weight(1,istsym),1)
     1    /dsum(nstate,weight(1,istsym),1)
      if (ipr.eq.2) vt = (pav-enext)/pav
160   continue
      do 170 istate=1,nstate
      if(weight(istate,istsym).ne.0d0) nstate_av=nstate_av+1
      dt(istate) = dsqrt(dt(istate))
      istat=iroot1+istate-1
      isy=istsy(istsym)
      isyref=isy
      write(iout,60) istat,isy
60    format(/' Results for state',i2,'.',i1/
     1        ' ====================='/)
      zekin=dmat(istate,2)
      en1=dmat(istate,3)
      en2=ener(istsym,istate)-en1
      en1=en1-potnuc
      vr=(zekin-ener(istsym,istate))/zekin
      call output_result('MCSCF','Energy',ener(istsym,istate)
     &   ,istate,isy)
      write (iout,61)
     1             potnuc,zekin,en1,en2,vr
61    format(
     1       1x,'Nuclear energy     ',t36,f15.8/
     1       1x,'Kinetic energy     ',t36,f15.8/
     1       1x,'One electron energy',t36,f15.8/
     1       1x,'Two electron energy',t36,f15.8/
     1       1x,'Virial ratio       ',t36,f15.8/)
      ityp=min0(ityp,2)
      call output_result('MCSCF','Dipole moment',dmat(istate,4:6)
     &   ,istate,isy,numberformat='3f15.8',debye=.TRUE.)
!      write (iout,1651) istate,isy,(dmat(istate,i),i=4,6),
!     1               (dmat(istate,i)*debye,i=4,6)
1651  format(' !MC STATE',i2,'.',i1,' DIPOLE MOMENTS:',t36,3f15.8,
     > ' a.u.',2x,3f12.6,' Debye')
      call setvar('ENERGY',ener(istsym,istate),'HARTREE',1,1,mxv,-1)
      call setvar('DMX',dmat(istate,4),'AU',1,1,mxv,-1)
      call setvar('DMY',dmat(istate,5),'AU',1,1,mxv,-1)
      call setvar('DMZ',dmat(istate,6),'AU',1,1,mxv,-1)
      do 165 i=1,n1elec
      ipr=i1elec(i)
      if(ipr.eq.1.or.ipr.gt.6) then
      write(iout,115) chars(ipr),dmat(istate,ipr)
115   format(1x,a20,t36,f15.8)
      end if
165   continue
170   continue
      call corlsr(iden)
200   continue
      if(many) then
        istate=-nstate_av
        name=dump_record
        ifil=dump_file
        iset=iset+1
        call write_den(q(idpos),iset,'CHARGE')
        write(iout,201) 'charge', name,ifil,iset
201     format(/' State-averaged ',a,' density matrix saved on record',
     1         i6,'.',i1,' (density set',i2,')')
        if(dmsav.gt.1.1d0) then
          if (dnrm2_X(ntdg,q(idpos_spin),1).gt.1d-3) then
c...   we managed to accumulate a spin density
            iset=iset+1
            call write_den(q(idpos_spin),iset,'SPIN')
            write(iout,201) 'spin', name,ifil,iset
          end if
        end if
        codxtr(3)=name*10+ifil
        codxtr(4)=0
        istate=istate_ref
      end if
      call corlsr(icivec)
      npr2=nact**2
      iden=igam+ic1d
      idipx=icorr(npr2)
      do 500 iop=1,noper1
      lexpec=iexp1(iop).ne.0
      str=oper1(iop)
      call nam_op(str,comp,ncomp,np,fak,center,ls)
      call op_center(iop,center,orig)
      do 490 i=1,ncomp
      first=.true.
      type=comp(i)
      ls=lenstr(type)
      if(ls.eq.0) goto 490
      if(type(1:2).eq.'DM') then
        if(nstsym.eq.1.and.nstate.eq.1) then
           goto 490
         else
           iofd=1
        end if
      end if
      call pgetq(q(idipx),type,znucx,zcorx,orig_op(1,iop))
      if(ls.gt.0) then
c.....clear variables for properties
        call setvarx(type(1:ls),0.0d0,' ',0)
        if(.not.lexpec) call setvarx('TR'//type(1:ls),0.0d0,' ',0)
      end if
      iofd=0
      istt=0
      do 350 istsym=1,nstsym
      isy=istsy(istsym)
      ivec1=icorr(ncisym(istsym)+lseg)
      do 300 istate=1,nstats(istsym)
      istt=istt+1
      call cget(q(ivec1),istsym,-istate)
      jst1=1
      if(lexpec) jst1=istsym
      do 295 jstsym=jst1,istsym
      jsy=istsy(jstsym)
      ivec2=icorr(ncisym(jstsym)+lseg)
      ja=1
      if(lexpec) ja=istate
      je=nstats(jstsym)
      if(istsym.eq.jstsym) je=istate-iofd
      do 290 jstate=ja,je
      call cget(q(ivec2),jstsym,-jstate)
      call dens1a (q(ivec1),q(ivec2),q(igam))
      trx= ddot_X(npr2,q(iden),1,q(idipx),1)*fak
      opnuc=0d0
      if(istate.eq.jstate.and.istsym.eq.jstsym) then
        trx=trx+zcorx*fak+znucx
        opnuc=znucx
      end if
      if(ls.gt.0) then
        if(istate.eq.jstate.and.istsym.eq.jstsym) then
          call setvarx(type(1:ls),trx,'AU',-1)
        else
          call setvarx('TR'//type(1:ls),trx,'AU',-1)
        end if
      end if
      if(dabs(trx).lt.1.d-8) goto 290
      if(first) write(iout,'()')
      first=.false.
      call prtran(trx,opnuc,istate,isy,jstate,jsy,type,'MCSCF',' ',orig)
290   continue
295   call corlsr(ivec2)
300   continue
350   call corlsr(ivec1)
c
490   continue
500   continue
      call corlsr (ibase)
      return
      end
      subroutine expec2
      implicit double precision (a-h,o-z)
      character*1 char1,char2
      character*16 type1,type2,type3,comp1,comp2,comp3
      character*16 center1,center2,center3
      logical expect
      include "common/tapes"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      include "common/clseg"
      include "common/cstate"
      include "common/ctran2"
      common /cjsym/ jstsym,isymd
      dimension expc(mxstt,mxstsy),orig1(3),orig2(3),orig3(3)
      call fzero(expc,mxstt*mxstsy)
c...  input of orig not supported
      call fzero(orig1,3)
      call fzero(orig2,3)
      call fzero(orig3,3)
      lxx=0
      lyy=0
      lzz=0
      nop=0
      do 1000 iop=1,ntran2
      if(nop.eq.0) write(iout,'()')
      nop = nop + 1
      expect=iexp2(iop).ne.0
      type1=tran2(1,iop)
      type2=tran2(2,iop)
      type3=tran2(3,iop)
      if(type1.ne.' ') call nam_op(type1,comp1,nc1,np1,fak1,center1,lc1)
      if(type2.ne.' ') call nam_op(type2,comp2,nc2,np2,fak2,center2,lc2)
      if(type3.ne.' ') call nam_op(type3,comp3,nc3,np3,fak3,center3,lc3)
      type1=comp1
      type2=comp2
      type3=comp3
      ls1=lenstr(type1)
      ls2=lenstr(type2)
      ls3=lenstr(type3)
c...  clear variables for properties
      call setvarx(type3(1:ls3),0.0d0,' ',0)
      call setvarx('TR'//type3(1:ls3),0.0d0,' ',0)
      if(type2.eq.' ') type2=type1
c...  get integrals
      n = nfreez+ncore+nact
      nc = nfreez+ncore
      n2 = n**2
      intgrl = icorr(n2)
      intgrr = icorr(n2)
      if(type3.ne.' ') intgr1 = icorr(nact**2)
      if(isign(1,np1).ne.isign(1,np2)) then
        write(iout,50) np1,np2
  50    format(' OPERATORS FOR EXPEC2 NOT OF EQUAL PARITY:',2i8)
        call fehler
      end if

      call pgetqc(q(intgrl),type1,znuc1,zcor1,orig1)
      call pgetqc(q(intgrr),type2,znuc2,zcor2,orig2)
      if(type3.ne.' ') then
        call pgetq (q(intgr1),type3,znuc3,zcor3,orig3)
      end if
c loop over states, constructing transition density
      igam = icorr(ne)
      igam1= igam+ic1d
      icivec=icorr(lenciv)
      do 350 istsym=1,nstsym
      jstsym=istsym
      nci=ncisym(istsym)
      nstate=nstats(istsym)
      call cget (q(icivec),istsym,nstate)
      do 300 istate=1,nstate
      ivec1=icivec+(istate-1)*nci
      ja=1
      if(expect) ja=istate
      do 290 jstate=ja,istate
      ivec2=icivec+(jstate-1)*nci
      call denst2 (q(ivec1),q(ivec2),q(igam))
      value=expecv(q(intgrl),q(intgrr),q(igam+ic1d),q(igam),n,nc)
      value=value*fak1*fak2
      if(np1.le.0) value=-value
      if(np3.ne.0) then
        value1=ddot_X(nact**2,q(intgr1),1,q(igam1),1)*fak3
        if(istate.eq.jstate) value1=value1+zcor3*fak3
        value=value+value1
      end if
      isy=istsy(istsym)
      jsy=istsy(jstsym)
      if(istate.eq.jstate.and.np1.eq.np2) then
        if(type3.eq.'LXLX'.and.lxx.eq.0) then
          expc(istate,istsym)=expc(istate,istsym)+value
        else if(type3.eq.'LYLY'.and.lyy.eq.0) then
          expc(istate,istsym)=expc(istate,istsym)+value
        else if(type3.eq.'LZLZ'.and.lzz.eq.0) then
          expc(istate,istsym)=expc(istate,istsym)+value
        end if
      end if
      if(istate.eq.jstate.and.istsym.eq.jstsym) then
        call setvarx(type3(1:ls3),value,'AU',-1)
      else
        call setvarx('TR'//type3(1:ls3),value,'AU',-1)
      end if
      call prtran(value,0.d0,istate,isy,jstate,jsy,type3,'MCSCF',' ',
     &   orig1)
290   continue
300   continue
350   continue
      if(type3.eq.'LXLX') lxx=1
      if(type3.eq.'LYLY') lyy=1
      if(type3.eq.'LZLZ') lzz=1
      call corlsr(intgrl)
1000  continue
      if(lxx*lyy*lzz.ne.0) then
        comp1='L**2'
        call clearvar('LL')
        do 250 istsym=1,nstsym
        do 250 istate=1,nstats(istsym)
        isy=istsy(istsym)
        call prtran(expc(istate,istsym),0.d0,istate,isy,istate,isy,
     1     comp1,'MCSCF',' ',orig1)
        call setvar('LL',expc(istate,istsym),'AU',1,1,mxm,-1)
250   continue
      end if
      return
      end
      function expecv(op1,op2,den1,den2,n,nc)
      implicit double precision (a-h,o-z)
c.....evaluates expecv=den2(ij,kl)*op1(ij)*op2(kl)
      include "common/cstate"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      dimension op1(n,n),op2(n,n),den1(nact,nact),den2(1)
c...  active-active contributions
      value=0d0
c     call outsqr(op1,n,n,n,'oper 1')
c     call outsqr(op2,n,n,n,'oper 2')
      do 10 i=1,nact
      do 10 j=1,nact
      vij = op1(nc+j,nc+i)
      ij = (i-1)*nact+j
      do 10 k=1,nact
      do 10 l=1,nact
      vkl = op2(nc+l,nc+k)
      kl = (k-1)*nact+l
      ijj=max(ij,kl)
      ijkl=ijj*(ijj-1)/2+min(ij,kl)
10    value = value + den2(ijkl)*vij*vkl
c...  active-core contributions, coulomb
      if(nc.eq.0) goto 60
      sumcor = 0d0
      sumcrr = 0d0
      do 20 i=1,nc
      sumcor = sumcor + op1(i,i)
20    sumcrr = sumcrr + op2(i,i)
      sumcor = 2d0*sumcor
      sumcrr = 2d0*sumcrr
      val=0d0
      do 30 k=1,nact
      do 30 l=1,nact
30    val=val+den1(l,k)*(sumcor*op2(nc+l,nc+k)
     1                  +sumcrr*op1(nc+l,nc+k))
      value=value+val
c...  active-core contributions, exchange
      val=0d0
      do 40 k=1,nact
      do 40 l=1,nact
      dkl=2.0d0*den1(k,l)
      do 40 i=1,nc
40    val = val - dkl*op1(nc+k,i)*op2(i,nc+l)
      value=value+val
c...  core-core contributions
      ov=dsum(nact,den1,nact+1)/dble(istnel(istsym))
      val = sumcor*sumcrr*ov
      value = value + val
      val=0d0
      do 50 i=1,nc
      do 50 j=1,nc
50    val = val - 2.0d0*op1(i,j)*op2(j,i)
      val=val*ov
      value=value+val
60    expecv=value
      end
      subroutine pgetq(p,oper,znuc,zcor,orig)
c...  reads a operator oper and transforms to MO
c...  p over active only, no symmetry pack, with off-diag
c...  core contribution in zcor
      implicit double precision (a-h,o-z)
      logical lcore
      character*(*) oper
      character*16 string
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/ctran2"
      include "common/big"
      include "common/tapes"
      include "common/cbas"
      dimension p(1),ivec(8),orig(3)
      lcore = .false.
      ndim = nact
      goto 10

      entry pgetqc(p,oper,znuc,zcor,orig)
      lcore = .true.
      ndim = nfreez+ncore+nact
10    continue
      iex1=iexcom_status()
      call excom(2)
      ixx=icorr(ntqg)
      iscr=icorr(2*ntqg)
c
c...  read orbitals
      do 20 isym=1,nirr
20    ivec(isym) = icorr(nsymao(isym)**2)
      call qget (ivec)
c
c
c...  copy active or occupied blocks to p
      call fzero (p,ndim*ndim)
      zcor=0d0
c
      ioff=0
      ipmo=icori(ndim)-1
      ismo=icori(ndim)-1
      ii=0
      if(lcore) then
        do is=1,nirr
          do i=1,ifreez(is)
            ii=ii+1
            iq(ipmo+ii)=i
            iq(ismo+ii)=is
          end do
        end do
        do is=1,nirr
          do i=1,ncor(is)
            ii=ii+1
            iq(ipmo+ii)=ifreez(is)+i
            iq(ismo+ii)=is
          end do
        end do
      end if
      do is=1,nirr
        do i=1,nactt(is)
          ii=ii+1
          iq(ipmo+ii)=ifreez(is)+ncor(is)+i
          iq(ismo+ii)=is
        end do
      end do
c...  read operator
      do 100 isy=1,nirr
      string=oper
      call read_op(q(ixx),znuc,string,0,isy,0,orig,isyop)
      if(isyop.eq.0) goto 100
c...  transform operator to MO basis
      call tranop(q(ixx),q(ivec(1)),q(iscr), isy,0, 1,0)
      ij=0
      do jj=1,ndim
        j=iq(ipmo+jj)
        js=iq(ismo+jj)
        do ii=1,ndim
          is=iq(ismo+ii)
          ij=ij+1
          if(mult(is,js).eq.isy) then
            ijj=ixx+ntqs(isy,is)+(j-1)*ntb(is)+iq(ipmo+ii)-1
            p(ij)=q(ijj)
          end if
        end do
      end do
      if(.not.lcore.and.isy.eq.1) then
        do is=1,nirr
          nic=ifreez(is)+ncor(is)
          if(nic.gt.0) then
            ijj=ixx+ntqs(isy,is)
            zcor=zcor+2.d0*dsum(nic,q(ijj),ntb(is)+1)
          end if
        end do
      end if
100   continue
      call corlsr (ixx)
      call excom(iex1)
      return
      end
      subroutine pget (p,iprop,znuc,zcor)
c...  p over active only, symmetry packed
c     zcor = 2*tr@over cor+fzc@(p)
      implicit double precision (a-h,o-z)
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      dimension p(1)
      dimension ix(8),ivec(8)
c
      do 10 isym=1,nirr
10    ix(isym) = icorr((nsymao(isym)*(nsymao(isym)+1))/2)
      call get1 (q(ix(1)),iprop,znuc)
c
      do 20 isym=1,nirr
20    ivec(isym) = icorr(nsymao(isym)**2)
      call qget (ivec)
c
      zcor=0d0
      ipp = 1
      do 50 isym=1,nirr
      ipp1=ipp
      n =nsymao(isym)
      m = nactt(isym)
      mc = ncor(isym) + ifreez(isym)
      l= m + mc
      if(l.eq.0) goto 50
      iy = icorr(n*n)
      iz = icorr(n*l)
      call expdm(q(ix(isym)),q(iy),n,1)
      call mxma(q(ivec(isym)),n,1,q(iy),1,n,q(iz),1,l,l,n,n)
      call mxma(q(iz),1,l,q(ivec(isym)),1,n,q(iy),1,l,l,n,l)
      zcor=zcor+2.0d0*dsum(mc,q(iy),l+1)
      do 45 i=1,m
      ii=(i+mc-1)*l+mc+iy-1
      do 45 j=1,i
      p(ipp)=q(ii+j)
45    ipp=ipp+1
      call corlsr (iy)
50    continue
      call corlsr (ix(1))
      return
      end
      subroutine eprint
      implicit double precision (a-h,o-z)
      character*8 tchar
      include "common/molen"
      include "common/etig"
      include "common/zahl"
      include "common/tapes"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/cstate"
      include "common/syminf"
      include "common/jobopt"
      include "common/cgeom"
      include "common/czmat"
      include "common/cwsave"
      include "common/big"
      include "common/cmpp"
      common/cener/ ener(mxstsy,mxstt)
      ibase = icorr(0)
      icivec = icorr(lenciv)
c     call allbuf(num6,-1)
      if(iguess.eq.0) call caserr('NO CI VECTOR PROVIDED')
      iwrnr=1
      do 100 istsym=1,nstsym
      nci=ncisym(istsym)
      nstate=nstats(istsym)
      call cget(q(icivec),istsym,nstate)
      ivec = icivec
      izint = icorr(ne)
      idiag = izint
      call fzero(q(izint),ne)
      call ci0 (q(izint),q(idiag))
      isig = icorr(nci)
      do 10 istate=1,nstate
      call fzero (q(isig),nci)
      call ci1 (q(ivec),q(isig),q(idiag),q(izint))
      ener(istsym,istate) = ddot_X(nci,q(ivec),1,q(isig),1)+core
10    ivec=ivec+nci
      call corlsr(izint)
      if(grstat.eq.z0) then
        do 20 i=10-nstate,1,-1
        tchar=gesnam(i)(1:min(len(gesnam(i)),len(tchar)))
        gesnam(i+nstate)= tchar
20      energ(i+nstate)=energ(i)
        do 30 istate=1,nstate
        if (istate.le.10) then
          gesnam(min(nstate,10)+1-istate)=befehl(1)
          energ(min(nstate,10)+1-istate)=ener(istsym,istate)
        endif
30      continue
      end if
100   continue
      istate=int(grstat)
      if(istate.le.0) istate=1
      istsym=int((grstat-dble(istate))*10.1d0)
      if(istsym.le.0) istsym=1
      if(igrty.le.6) then
        zenergy=ener(istsym,istate)
        if(igrty.eq.6) then
          istate=int(grstat2)
          if(istate.le.0) istate=1
          istsym=int((grstat2-dble(istate))*10.1d0)
          if(istsym.le.0) istsym=1
          zenergy=ener(istsym,istate)-zenergy
          zmethod='DEMC'
        end if
        if(grstat.ne.z0) then
          do i=9,1,-1
            gesnam(i+1)=gesnam(i)
            energ(i+1)=energ(i)
          end do
          gesnam(1)=befehl(1)
          energ(1)=zenergy
        end if
      else if(igrty.gt.6) then
        zenergy=0d0
        zmethod='NACME'
      end if
      call geomw(0)
      call corlsr(ibase)
      return
      end
      subroutine ci0 (zint,diag)
      implicit double precision(a-h,o-z)
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      integer ijpos(8),ikpos(8)
      dimension zint(ne),diag(nci)
c statement function ind is not used in the current subroutine.
c      ind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      call citime('ci0',0)
c
c...  load integrals
      ibase = icorr(0)
      do 10 isym=1,nirr
      ijpos(isym) = icorr(maxbas**2)
10    ikpos(isym) = icorr(maxbas**2)
      do 20 i=1,nprim
      do 20 j=1,i
20    call loadjk (ijpos,ikpos,i,j,zint)
      call loadjk (ijpos,ikpos,0,0,zint)
      call corlsr (ibase)
c
 999  call citime('ci0',1)
      return
      end
      subroutine print(iprin)
c
c     printing routine. called by anal (lprint) ? wvfn (iprint)
c
      implicit double precision(a-h,o-z)
      include "common/tapes"
      include "common/maxbfn"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      include "common/cref"
      include "common/cstate"
      logical btest
      include "common/detcic"
      character*9 word
      dimension iqpos(8),occno(31),scra(31),orbnat(31)
c
      ibase = icorr(0)
c
      ifzcpr=0
      if (btest(iprin,12).or.btest(iprin,10)) then
      icivec=icorr(lenciv)
      if(iguess.eq.0) call caserr('NO CI VECTOR PROVIDED')
      end if
      if (btest(iprin,9).or.btest(iprin,12).or.btest(iprin,13)) then
      do 20 isym=1,nirr
20    iqpos(isym) = icorr(nsymao(isym)**2)
      call qget (iqpos)
      end if
c
c...  =orbitals=
      if (btest(iprin,9)) then
      write(iout,140)
140   format(/' MOLECULAR ORBITAL COEFFICIENTS'
     1       /' =============================='/
     2        '  NR SYM      TYP       COEFFICIENTS')
      do 170 isym=1,nirr
      n=min0(nsymao(isym),ifreez(isym)+nprm(isym)+2)
      ia=iqpos(isym)-1
      do 170 i=1,n
      if (ifzcpr.ne.0.and.i.le.ifreez(isym)) goto 170
      word='VIRTUAL'
      if (i.le.ifreez(isym)+nprm(isym))  word ='ACTIVE '
      if (i.le.ifreez(isym)+ncor(isym))  word ='CORE   '
      if (i.le.ifreez(isym))             word ='FROZEN '
      write(iout,160) i,isym,word,(q(ia+k),k=1,nsymao(isym))
170   ia=ia+nsymao(isym)
160   format(/1x,2i3,5x,a7,(t20,10f11.6))
      ifzcpr=1
      end if
c
c... =natorb=
      if (btest(iprin,12)) then
      igam = icorr(ne)
      do 270 istsym=1,nstsym
      nci=ncisym(istsym)
      nstate=nstats(istsym)
      call cget(q(icivec),istsym,nstate)
      ivec=icivec
      if(nstsym.gt.1) write(iout,201) istsym
201   format(/' STATE SYMMETRY',i2)
      do 260 istat1=1,nstate
      do 260 istat2=1,istat1
      call denst1(q(ivec+(istat1-1)*nci),q(ivec+(istat2-1)*nci),q(igam))
      if (istat1.eq.istat2) then
      write(iout,180) istat1
180   format(/' ONE-PARTICLE DENSITY MATRIX FOR STATE',i2/
     1      ' =======================================')
      call outsqr (q(igam+ic1d),nact,nact,nact,' ')
c...  natural orbital analysis
      write(iout,190)
190   format(/' NATURAL ORBITALS'/' ================')
      call f02abf (q(igam+ic1d),nact,nact,occno,q(igam),nact,scra,0)
      do 240 i=1,nact
      occ = occno(nact+1-i)
      call fmove (q(igam+(nact-i)*nact),orbnat,nact)
      zmax = 0.0d0
      do 200 j=1,nact
      if (dabs(orbnat(j)).le.zmax) goto 200
      max = j
      zmax = dabs(orbnat(j))
200   continue
      isym = itypea(max)
      ig = icorr(nsymao(isym))
      if (zmax.lt.0.0d0) call dscal_X (nact,-1.0d0,orbnat,1)
      kk=ncor(isym)+ifreez(isym)
      call fzero (q(ig),nsymao(isym))
      do 210 j=1,nact
      if (itypea(j).ne.isym) goto 210
      call daxpy_X(nsymao(isym),orbnat(j),
     1             q(iqpos(isym)+kk*nsymao(isym)),1,q(ig),1)
      kk = kk + 1
210   continue
      write(iout,220) i+nfreez+ncore,occ,isym,(orbnat(j),j=1,nact)
      write(iout,230) (q(ig-1+j),j=1,nsymao(isym))
220   format(/i4,f10.7,1x,'SYM',i2,2x,
     1'E-VEC',8(t29,9f11.7/) )
230   format(/22x,'COEFFS',70(t29,9f11.7/) )
      call corlsr(ig)
240   continue
      else
      write(iout,250) istat1,istat2
250   format(/' ONE-PARTICLE TRANSITION DENSITY MATRIX FOR STATES',2i2/
     1      ' =====================================================')
      call outsqr (q(igam+ic1d),nact,nact,nact,' ')
      end if
260   continue
270   continue
      call corlsr (igam)
      end if
c
c... =civector=
      if (btest(iprin,10)) then
      call clearvar('COEFF')
      do 300 istsym=1,nstsym
      nci=ncisym(istsym)
      nstate=nstats(istsym)
      call cget(q(icivec),istsym,nstate)
      call prici(icivec)
300   continue
      end if
c
c... =formulae=
      if (btest(iprin,2)) then
      end if
c
      call corlsr (ibase)
c
      return
      end
      subroutine savdm(gam,denav,orb,iset,type,save)
      implicit double precision(a-h,o-z)
      logical save
      include "common/maxatm"
      include "common/maxbfn"
      include "common/dumpinfow"
      include "common/code"
      include "common/cstate"
      include "common/corbdim"
      include "common/casscf"
      include "common/syminf"
      include "common/jobopt"
      include "common/big"
      include "common/cref"
      include "common/cbas"
      include "common/tapes"
      dimension gam(nact,nact),orb(1),denav(*)
      character *(*) type
      character *16 typ
      call excom(2)
      npri = nprim+nfreez
      igam=icorr(npri**2)
      iden=icorr(ntqg)
      iscr=icorr(max(ntdgc(1),ntqg))
      do 40 isym=1,nirr
      nc=ifreez(isym)+ncor(isym)
      npri=nc+nactt(isym)
      idpos=iden+ntq(isym)
      if(npri.eq.0) then
        call fzero(q(idpos),nt(isym)**2)
        goto 40
      end if
      call fzero(q(igam),npri**2)
      ii=igam
      typ=type
      call upcase(typ)
      if (typ.eq.'CHARGE') then
       occ_core=2d0
      else if (typ.eq.'SPIN') then
       occ_core=0d0
      else
       write (6,*) typ
       call Error('Unknown density type','muprop:savdm')
      end if
      do 10 i=1,nc
      q(ii)=occ_core
10    ii=ii+npri+1
      ii=nc
      do 30 i=1,nact
      if (itypea(i).ne.isym) goto 30
      ii=ii+1
      jj=nc
      do 20 j=1,nact
      if (itypea(j).ne.isym) goto 20
      jj=jj+1
      q(igam+(jj-1)+(ii-1)*npri) = gam(j,i)
20    continue
30    continue
      call mxma(orb(ntq(isym)+1),1,ntb(isym),q(igam),1,npri,
     1          q(iscr),1,ntb(isym),nt(isym),npri,npri)
      call mxma(q(iscr),1,ntb(isym),orb(ntq(isym)+1),ntb(isym),1,
     1          q(idpos),1,ntb(isym),nt(isym),npri,nt(isym))
40    continue
      call reduc(q(iden),q(iscr),1,1,1.0d0,0)
      call daxpy_X(ntdg,weight(istate,istsym),q(iscr),1,denav,1)
      if(save) then
        name=dump_record
        ifil=dump_file
        isyref_save=isyref
        isyref=istsy(istsym)
        iset=iset+1
        call write_den(q(iscr),iset,typ)
        isyref=isyref_save
        call locase(typ)
        write(iout,50) typ(1:lenstr(typ)),istate,istsy(istsym),
     1                 name,ifil,iset
50      format(' First order ',a,' density matrix for state',i2,'.',i1,
     1     ' saved on record',i6,'.',i1,' (density set',i2,')')
        codxtr(3)=name*10+ifil
        codxtr(4)=0
      end if
      call excom(1)
      call corlsr(igam)
      return
      end
