package HackaMol::X::NERF;
 # ABSTRACT: Natural extension reference frame implementation for molecular building 

use 5.008;
use Moose;
use namespace::autoclean;
use Carp;
use Math::Vector::Real;
use Math::Trig; 

with 'HackaMol::Roles::NERFRole';

has 'bb_dihes' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[ArrayRef]',
    default => sub { [] },
    lazy    => 1,
    handles => {
      has_bb_dihes => 'count',
      get_bb_dihes => 'get',
      set_bb_dihes => 'set',
      all_bb_dihes => 'elements',
    }, #[$phi,psi,omega] N-CA---C-N; CA-C---N-CA; C-N---CA-C
    # thus there are only bb 
);


__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

       #Let's build a six member ring 
       use HackaMol::X::NERF;
    
       my $bld = HackaMol::X::NERF->new;
       
       push @vecs, $bld->init() ; 
       push @vecs, $bld->extend_a(  $vecs[0]  ,   1.47              );
       push @vecs, $bld->extend_ab( @vecs[0,1],   1.47, 109.5       );
       push @vecs, $bld->extend_abc(@vecs[0,1,2], 1.47, 109.5,  60 );
       push @vecs, $bld->extend_abc(@vecs[1,2,3], 1.47, 109.5, -60 );
       push @vecs, $bld->extend_abc(@vecs[2,3,4], 1.47, 109.5,  60 );

       printf ("C %10.6f %10.6f %10.6f\n", @$_ ) foreach @vecs;

=head1 DESCRIPTION

The HackaMol::X::NERF library (HMX::NERF) consumes HackaMol::Roles::NERFRole, which is provided with the 
HackaMol core.  HM::Roles::NERFRole provides a quick implementation of the Natural Extension Reference Frame 
method for building cartesian coordinates from internal coordinates.  HMX::NERF seeks to provide methods 
to make building molecules and peptides more convenient.  

The API is changing/expanding.  

=method build_extended_peptide ('aaaaaaaaapaaaaapaaaaa')

take a string of amino acids (1-letter code) and build an exteded backbone and sidechains. 
phi, psi, and omega: -120,140,180  
Prolines are added wih -60, 140,180.


=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<HackaMol::AtomGroup>
* L<HackaMol::Molecule>
* L<HackaMol::Roles::NERFRole>
 
