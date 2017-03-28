MODULE pdf_calculator

    USE double
    IMPLICIT NONE

    public :: calc_pdf

CONTAINS

SUBROUTINE calc_pdf(positions, cell, num_atoms, pdf, r_space, dr, num_images)

    !--------------------------------------------------------------------------!
    ! Purpose:                                                                 !
    ! This routine takes MD data as input and produces a PDF.                  !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! positions(in)   : Array of MD trajectories of shape (3, n_atoms, n_t)    !
    ! cell(in)        : Array of lattice vectors of shape (3, 3, n_t)          !
    ! pdf(inout)      : pdf normalised with 1/(4*PI*r^2*num_pairs*dr)          !
    ! num_atoms(in)   : Number of atoms in unit cell                           ! 
    ! r_space(inout)  : Array of distances for which to calculate pdf          !
    ! dr (inout)      : Resolution or bin width of pdf in Angstrom             !
    !--------------------------------------------------------------------------!
    ! Modules used: double                                                     !
    !--------------------------------------------------------------------------!

    ! radial parameters
    real(kind=dp), intent(in) :: dr
    integer                   :: r_length
    ! array sizes
    integer, intent(in) :: num_atoms
    ! data arrays
    real(kind=dp), intent(in), dimension(:,:)  :: positions
    real(kind=dp), intent(in), dimension(3,3)  :: cell
    real(kind=dp), intent(inout), dimension(:)   :: pdf
    real(kind=dp), intent(inout), dimension(:)   :: r_space
    real(kind=dp), dimension(3) :: translate
    ! translation vector magnitude
    real(kind=dp)  :: trans_mag     = 0.0
    ! number of images from each edge
    integer, intent(inout) :: num_images
    ! error codes, indices and utils
    logical :: debug    = .false.
    integer :: iion, jion, kpos = 0
    integer :: a, b, c  ! coefficents of LV's for periodic images
    integer :: bindex
    real(kind=dp)  :: tmp_dr

    ! A L L O C A T I O N

    pdf = 0.0_dp
    r_length = size(r_space)
    do iion=1, r_length
        r_space(iion) = iion*dr
    end do

    ! C A L C U L A T I O N

    do iion=1,num_atoms
        do jion=iion+1,num_atoms
            tmp_dr = 0.0_dp
            do kpos=1,3
                tmp_dr = tmp_dr + (positions(kpos, iion) - positions(kpos, jion))**2
            end do
            tmp_dr = sqrt(tmp_dr)
            bindex = nint(tmp_dr/dr)
            if(bindex>r_length) cycle
            ! increment pdf by 2 as zeroth cell has double weight of images
            pdf(bindex) = pdf(bindex) + 2/(num_atoms * (num_images+1)**3)
            if(debug) then
                if(iion.eq.1) write(*,*) tmp_dr, ' going to the bin for ', r_space(bindex), '->', r_space(bindex+1)
            end if
        end do
    end do

    ! P E R I O D I C

    if(num_images>0) then
        do a=-1*num_images, num_images
            do b=-1*num_images, num_images
                do c=-1*num_images, num_images
                    trans_mag = 0.0_dp
                    do kpos=1,3
                        translate(kpos) = a * cell(kpos, 1)
                        translate(kpos) = translate(kpos) + b * cell(kpos, 2)
                        translate(kpos) = translate(kpos) + c * cell(kpos, 3)
                    end do
                    if(a.eq.0 .and. b.eq.0 .and. c.eq.0) cycle
                    do iion=1,num_atoms
                        do jion=1,num_atoms
                            tmp_dr = 0.0_dp
                            do kpos=1,3
                                tmp_dr = tmp_dr + (positions(kpos, iion) - positions(kpos, jion) + translate(kpos))**2
                            end do
                            tmp_dr = sqrt(tmp_dr)
                            bindex = nint(tmp_dr/dr)
                            if(bindex>r_length) cycle
                            pdf(bindex) = pdf(bindex) + 1/(num_atoms * (num_images+1)**3)
                        end do
                    end do
                end do
            end do
        end do
    end if
    write(*,*) 'pdf complete!'
    return

END SUBROUTINE calc_pdf

END MODULE pdf_calculator
