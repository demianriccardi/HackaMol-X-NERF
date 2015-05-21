#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;
use HackaMol::X::NERF;
use Math::Vector::Real;
use Data::Dumper;

#create an extension builder
my $bld = HackaMol::X::NERF->new;

#print Dumper $bld; exit;

my $mol = HackaMol->new->read_file_mol(shift);

my $new_mol = hbuild_sp3($mol,4);

sub hbuild_sp3{
  my ($mol, $iat) = (shift,shift);
  my $atom = $mol->get_atoms($iat);
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

