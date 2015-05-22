#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;
use HackaMol::X::NERF;
use Math::Vector::Real;
use Data::Dumper;

#create an extension builder
my $bld  = HackaMol::X::NERF->new;
my $hack = HackaMol->new;

#print Dumper $bld; exit;

my $mol = HackaMol->new->read_file_mol(shift);

my $new_mol = hbuild_sp3($mol,4);

sub hbuild_sp3{

  my ($mol, $iat) = (shift,shift);

  #get atom of interest
  my $atom = $mol->get_atoms($iat);
 
  #find all the bonds to the atom of interest 
  my @bonds = $hack->find_bonds_brute(
                                  bond_atoms => [$atom],
                                  candidates => [$mol->all_atoms],
                                  fudge      => 0.45,
                                  max_bonds  => 6,
  );
 
  $mol->push_bonds(@bonds);
  #my $number_of_hyd = 4 - $atom->bond_count;

  #add some number of hydrogens depending on number of bonds to atom of interest.
  if ($atom->bond_count == 1){
    my $at_bound = $bonds[0]->get_atoms(1);

    #find all bonds from the bound atom to establish the reference frame
    my @b_bonds = $hack->find_bonds_brute(
                                  bond_atoms => [$at_bound],
                                  candidates => [$mol->all_atoms],
                                  fudge      => 0.45,
                                  max_bonds  => 6,
    );
    # this may be cheating.
    pop @b_bonds;

    # map out the math::vector::real xyz for each of the atoms of interest 
    my ($a,$b,$c) = map{$_->get_atoms(1)->xyz} (@bonds,@b_bonds);
    
    # return a math::vector::real xyz for each of the 3 hydrogens
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109, 0);
    my $h2 = $bld->extend_abc($h1,$a, $atom->xyz, 1.09, 109, 120);
    my $h3 = $bld->extend_abc($h1,$a, $atom->xyz, 1.09, 109, -120);

    $mol->push_atoms(
                      map{HackaMol::Atom->new(Z => 1, coords => [$_])} ($h1,$h2,$h3)
    );  
    $mol->print_xyz;
  }

  exit;
 
  print Dumper $mol;exit;
  #print Dumper V(1,1,0) -> versor; exit;
  my $h1 = HackaMol::Atom->new(
                                symbol => 'H',
                                coords => [
                     $bld->extend_a($atom->xyz,'1.09', V(1,1,1))
                                ] 
                               );  
  my $new_mol = HackaMol::Molecule->new(atoms=>[$atom, $h1]);
  $new_mol->print_xyz;
  say $h1->distance($atom);

#print Dumper $h1;
#print Dumper $atom;
  
}

