module blocktypes

type Y01BlockData
     integer :: n
     integer :: l1, l2
     double precision,allocatable :: vec0(:)
end type Y01BlockData

type EblockData
integer :: n
integer :: l1,l2
integer,allocatable :: pos(:)
integer,allocatable :: ipiv(:)
double precision,allocatable :: vec(:), matX(:,:),matY(:,:)
end type

end module blocktypes
