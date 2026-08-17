module geom_input
      !
      ! Building the molecular system and the AO basis from the GammCor
      ! input file. All modules which need TSystem/TAOBasis before the
      ! integrals are computed should go through geom_ReadSystemBasis
      ! instead of calling sys_Read_XYZ and basis_init themselves.
      !
      use gammcor_integrals
      use types

      implicit none

      private
      public :: geom_ReadSystemBasis
      public :: xyz_source
      public :: geom_DeleteScratch

      character(*), parameter :: GEOM_INPUT_FILE = "./input.inp"
      character(*), parameter :: GEOM_SCRATCH_FILE = "gammcor_geom_tmp.inp"

contains

      subroutine geom_ReadSystemBasis(System, AOBasis, Flags, SortAngularMomenta, Units, InputFile)
            !
            ! Read the geometry pointed at by the GammCor input and set up
            ! the AO basis for it.
            !
            type(TSystem), intent(out)         :: System
            type(TAOBasis), intent(out)        :: AOBasis
            type(FlagsData), intent(in)        :: Flags
            logical, intent(in) :: SortAngularMomenta
            integer, optional, intent(in)      :: Units
            character(*), optional, intent(in) :: InputFile


            character(:), allocatable :: XYZPath
            integer :: Units0
            integer :: ShellOrder

            if (present(Units)) then
                  Units0 = Units
            else
                  Units0 = SYS_UNITS_ANGSTROM
            end if
            if (present(InputFile)) then
                  call xyz_source(XYZPath, InputFile)
            else
                  call xyz_source(XYZPath, GEOM_INPUT_FILE)
            end if
            call msg("Geometry read from " // XYZPath)

            call auto2e_init()
            call sys_Read_XYZ(System, XYZPath, Units0)

            if (SortAngularMomenta) then
               ShellOrder = SHELL_ORDER_BY_MOMENTUM
            else
               ShellOrder = SHELL_ORDER_FIXED
            end if

            if (Flags%BasisAssign%Initialized) then
                  call warn_basis_keyword_ignored(Flags)
                  call basis_init(AOBasis, System, &
                        ShellOrder=ShellOrder, &
                        BasisAssign=Flags%BasisAssign)
            else
                  call basis_init(AOBasis, System, &
                        FilePath=Flags%BasisSetPath // Flags%BasisSet, &
                        ShellOrder=ShellOrder)
            end if
            call AObasis%display()
      end subroutine geom_ReadSystemBasis

      subroutine warn_basis_keyword_ignored(Flags)
            !
            ! Say out loud that the BasisAssignement block takes precedence over
            ! the Basis keyword left in the Calculation block. Printed next to the
            ! basis set summary, so that the basis actually used in the calculation
            ! cannot be confused with the one given by the Basis keyword.
!
            use display, only: msg, MSG_WARNING
            implicit none
            type(FlagsData), intent(in) :: Flags
            
            if(.not.Flags%BasisAssign%Initialized) return
            if(.not.allocated(Flags%BasisSet)) return
            if(len_trim(Flags%BasisSet)==0) return
            
            call msg("Warning: the BasisAssignement block overrides the Basis keyword &
                  &from the Calculation block: " // trim(Flags%BasisSet) // " is ignored", &
                  MSG_WARNING)
            
      end subroutine warn_basis_keyword_ignored

      subroutine xyz_source(XYZFile, InputFile)
            !
            ! Decide which file sys_Read_XYZ should be pointed at.
            !
            ! The geometry can be given inline in the GammCor input,
            !
            !       xyz
            !       2
            !       N  0.0  0.0  0.0
            !       N  0.0  0.0  2.075
            !       end
            !
            ! or the xyz directive can point to an external file,
            !
            !       xyz my_geom.xyz
            !       xyz /home/.../my_geom.xyz
            !
            ! The external file may either repeat the "xyz ... end" block
            ! syntax (then it is passed to sys_Read_XYZ as it is) or be
            ! a plain, standard xyz file (NAtoms / comment / atom lines).
            ! In the latter case the file is rewritten into the block
            ! syntax in a scratch file, whose name is returned.
            !
            character(:), allocatable, intent(out) :: XYZFile
            character(*), intent(in)               :: InputFile

            character(:), allocatable :: ExtFile

            call xyz_ReadDirective(ExtFile, InputFile)
            if (len(ExtFile) == 0) then
                  !
                  ! No file name after the xyz keyword: the coordinates
                  ! are in the input file itself
                  !
                  XYZFile = InputFile
            else if (xyz_HasBlock(ExtFile)) then
                  XYZFile = ExtFile
            else
                  XYZFile = GEOM_SCRATCH_FILE
                  call xyz_ConvertPlain(ExtFile, XYZFile)
            end if
      end subroutine xyz_source

      subroutine geom_DeleteScratch()
            !
            ! Remove the scratch file written by xyz_ConvertPlain. Called
            ! at the end of the run: the same scratch file is rewritten and
            ! reread by every xyz_source call, so it cannot be deleted any
            ! earlier. Nothing happens if no scratch file has been written.
            !
            integer :: u
            logical :: FileExists

            inquire(file=GEOM_SCRATCH_FILE, exist=FileExists)
            if (FileExists) then
                  open(newunit=u, file=GEOM_SCRATCH_FILE, status="OLD")
                  close(u, status="DELETE")
            end if
      end subroutine geom_DeleteScratch

      subroutine xyz_ReadDirective(ExtFile, InputFile)
            !
            ! Extract the file name following the xyz keyword. An empty
            ! string is returned if the xyz keyword is absent or is not
            ! followed by a file name.
            !
            character(:), allocatable, intent(out) :: ExtFile
            character(*), intent(in)               :: InputFile

            integer :: u
            logical :: eof
            character(:), allocatable :: line, key, val

            ExtFile = ""
            u = io_text_open(InputFile, "OLD")
            lines: do
                  call io_text_readline(line, u, eof)
                  if (eof) then
                        exit lines
                  end if
                  call split(line, key, val)
                  if (uppercase(key) == "XYZ") then
                        !
                        ! A bare xyz keyword opens an inline geometry block
                        ! and leaves ExtFile empty
                        !
                        if (.not. isblank(val)) then
                              ExtFile = trim(adjustl(val))
                        end if
                        exit lines
                  end if
            end do lines
            close(u)
      end subroutine xyz_ReadDirective

      function xyz_HasBlock(FilePath)
            !
            ! Check if FilePath uses the "xyz ... end" block syntax
            ! understood by sys_Read_XYZ.
            !
            logical :: xyz_HasBlock
            character(*), intent(in) :: FilePath

            integer :: u
            logical :: eof
            character(:), allocatable :: line, key, val

            xyz_HasBlock = .false.
            u = io_text_open(FilePath, "OLD")
            lines: do
                  call io_text_readline(line, u, eof)
                  if (eof) then
                        exit lines
                  end if
                  if (isblank(line) .or. iscomment(line)) then
                        cycle lines
                  end if
                  call split(line, key, val)
                  if (uppercase(key) == "XYZ") then
                        xyz_HasBlock = .true.
                        exit lines
                  end if
            end do lines
            close(u)
      end function xyz_HasBlock

      subroutine xyz_ConvertPlain(PlainFile, BlockFile)
            !
            ! Rewrite a standard xyz file
            !
            !       NAtoms
            !       comment line
            !       El  x  y  z
            !       ...
            !
            ! into the "xyz ... end" block expected by sys_Read_XYZ: the
            ! keywords are added and the comment line is dropped. Only a
            ! single molecule can be defined this way; the charge and the
            ! multiplicity keep their default values.
            !
            character(*), intent(in) :: PlainFile
            character(*), intent(in) :: BlockFile

            integer :: u, v, k, NAtoms, NRead, ios
            character(len=256) :: line
            logical :: eof

            open(newunit=u, file=PlainFile, status="OLD", action="READ")
            read(u, *, iostat=ios) NAtoms
            if (ios /= 0 .or. NAtoms <= 0) then
                  call msg("Could not read the number of atoms from the xyz file " &
                        // PlainFile, MSG_ERROR)
                  error stop
            end if
            open(newunit=v, file=BlockFile, status="REPLACE", action="WRITE")
            write(v, "(A)") "xyz"
            write(v, "(I0)") NAtoms
            NRead = 0
            call xyz_ReadNonBlank(u, line, eof)
            if (.not. eof) then
                  if (xyz_IsAtomLine(line)) then
                        write(v, "(A)") trim(line)
                        NRead = 1
                  end if
            end if
            do k = NRead + 1, NAtoms
                  call xyz_ReadNonBlank(u, line, eof)
                  if (eof) then
                        call msg("The xyz file " // PlainFile // " declares more atoms &
                              &than it contains", MSG_ERROR)
                        error stop
                  end if
                  write(v, "(A)") trim(line)
            end do
            write(v, "(A)") "end"
            close(v)
            close(u)
      end subroutine xyz_ConvertPlain

      subroutine xyz_ReadNonBlank(u, line, eof)
            !
            ! Read the next non-blank line of an open text file.
            !
            integer, intent(in)        :: u
            character(*), intent(out)  :: line
            logical, intent(out)       :: eof

            integer :: ios

            do
                  read(u, "(A)", iostat=ios) line
                  if (ios /= 0) then
                        line = ""
                        eof = .true.
                        return
                  end if
                  if (len_trim(line) > 0) then
                        eof = .false.
                        return
                  end if
            end do
      end subroutine xyz_ReadNonBlank

      function xyz_IsAtomLine(line)
            !
            ! Check if line is "El x y z", i.e. a known element symbol
            ! followed by three numbers, as opposed to a comment line.
            !
            logical :: xyz_IsAtomLine
            character(*), intent(in) :: line

            character(len=64) :: Element
            real(F64) :: R(3)
            integer :: ios

            xyz_IsAtomLine = .false.
            read(line, *, iostat=ios) Element, R(1), R(2), R(3)
            if (ios /= 0) then
                  return
            end if
            xyz_IsAtomLine = (znumber_short(trim(Element)) > 0)
      end function xyz_IsAtomLine
end module geom_input
