#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;
use HackaMol::X::NERF;
use Math::Vector::Real;
use Carp;
use Data::Dumper;

#create an extension builder
my $bld  = HackaMol::X::NERF->new;
my $hack = HackaMol->new;

#print Dumper $bld; exit;

#my $mol = HackaMol->new->read_file_mol('zmatrices/val.zmat');
my $mol = HackaMol->new->pdbid_mol('2cba');
#hbuild_sp3($mol,4);
#hbuild_sp3($mol,5);
#hbuild_sp3($mol,6);

my @fun = $mol->select_atoms( sub{
                                    $_->resname eq 'VAL'
                                   # or 
                                   # $_->resname eq 'ALA'
                                 });

my @groups =  $hack->group_by_atom_attr('resid',@fun);

my $newmol = HackaMol::Molecule->new(groups=>[@groups]);

#foreach my $iat (1 .. $mol->count_atoms){
#  hbuild_sp3($mol,$iat-1);
#}

foreach my $iat ( 1 .. $newmol->count_atoms ){
  hbuild_sp2($newmol,$iat-1);
}

my @hyd = grep { $_->Z == 1 } $newmol->all_atoms;
$mol->push_atoms(@hyd);
$mol->print_pdb('shit.pdb');

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

  say scalar(@bonds); say $mol->count_atoms;

  if (@bonds > 3){
    print "already at least four bonds on atom $iat\n";
    return;
  }
  elsif(@bonds == 1){ 
  #add some number of hydrogens depending on number of bonds to atom of interest.
    my $at_bound = $bonds[0]->get_atoms(1);

    #find all bonds from the bound atom to establish the reference frame
    my @b_bonds = $hack->find_bonds_brute(
                                  bond_atoms => [ $at_bound ],
                                  candidates => [ $mol->select_atoms( sub{ $_->iatom != $atom->iatom } ) ],
                                  fudge      => 0.45,
                                  max_bonds  => 6,
    );

    # map out the math::vector::real xyz for each of the atoms of interest 
    my ($a,$b,$c) = map{$_->get_atoms(1)->xyz} (@bonds,@b_bonds);
    print Dumper $a, $b,$c; 
   # return a math::vector::real xyz for each of the 3 hydrogens
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109, 0);
    my $h2 = $bld->extend_abc($h1,$a, $atom->xyz, 1.09, 109, 120);
    my $h3 = $bld->extend_abc($h1,$a, $atom->xyz, 1.09, 109, -120);

    $mol->push_atoms(
                      map{HackaMol::Atom->new(name=> 'HH', Z => 1, coords => [$_])} ($h1,$h2,$h3)
    );
    return;  
  }
  elsif(@bonds == 2){
    my ($a,$b) = map {$_->get_atoms(1)->xyz} @bonds;
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109,  120);
    my $h2 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109, -120);
    $mol->push_atoms(
                      map{HackaMol::Atom->new(name=> 'HH',Z => 1, coords => [$_])} ($h1,$h2)
    );
    return $mol;
  }
  elsif(@bonds == 3){
    my ($a,$b,$c) = map {$_->get_atoms(1)->xyz} @bonds;
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109,  120);
    if ($h1->dist($c) < 1.){
      $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109,  -120);
    }
    $mol->push_atoms( HackaMol::Atom->new(name=> 'HH', Z => 1, coords => [$h1]) );
    return $mol;
  }
  
}
# example valine
#hbuild_isp3_iaibic($mol,4,1,5,6)

sub hbuild_isp3_iaibic{

  my ($mol, $iat, $ia,$ib,$ic) = @_ ;

  croak "pass at least two indices of bound atoms" unless defined($ib);

  #get atom of interest
  my $atom = $mol->get_atoms($iat);
  my $ata = $mol->get_atoms($ia);
 
  if (defined ($ic)){
    my $atb = $mol->get_atoms($ib);
    my $atc = $mol->get_atoms($ic);
 
    my $h1 = $bld->extend_abc($atb->xyz, $ata->xyz, $atom->xyz, 1.09, 109,  120);
    if ($h1->dist($atc->xyz) < 1.){
      $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 109,  -120);
    }
    $mol->push_atoms( HackaMol::Atom->new(Z => 1, coords => [$h1]) );
    return $mol;
  }
  elsif (defined ($ib)){
    my $atb = $mol->get_atoms($ib);
    my $h1 = $bld->extend_abc($atb->xyz, $ata->xyz, $atom->xyz, 1.09, 109,   120);
    my $h2 = $bld->extend_abc($atb->xyz, $ata->xyz, $atom->xyz, 1.09, 109,  -120);
    $mol->push_atoms( HackaMol::Atom->new(Z => 1, coords => [$h1,$h2]) );
    return $mol;
  }
  
}

sub hbuild_sp2{

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

  say scalar(@bonds); say $mol->count_atoms;

  if (@bonds > 2){
    print "already at least three bonds on atom $iat\n";
    return;
  }
  elsif(@bonds == 1){ 
  #add some number of hydrogens depending on number of bonds to atom of interest.
    my $at_bound = $bonds[0]->get_atoms(1);

    #find all bonds from the bound atom to establish the reference frame
    my @b_bonds = $hack->find_bonds_brute(
                                  bond_atoms => [ $at_bound ],
                                  candidates => [ $mol->select_atoms( sub{ $_->iatom != $atom->iatom } ) ],
                                  fudge      => 0.45,
                                  max_bonds  => 6,
    );

    # map out the math::vector::real xyz for each of the atoms of interest 
    my ($a,$b) = map{$_->get_atoms(1)->xyz} (@bonds,@b_bonds);
   # return a math::vector::real xyz for each of the 3 hydrogens
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 120, 180);
    my $h2 = $bld->extend_abc($h1,$a, $atom->xyz, 1.09, 120, 180);

    $mol->push_atoms(
                      map{HackaMol::Atom->new(name=> 'HH', Z => 1, coords => [$_])} ($h1,$h2)
    );
    return;  
  }
  elsif(@bonds == 2){
    my ($a,$b) = map {$_->get_atoms(1)->xyz} @bonds;
    my $h1 = $bld->extend_abc($b, $a, $atom->xyz, 1.09, 120,  180);
    $mol->push_atoms(
                      map{HackaMol::Atom->new(name=> 'HH',Z => 1, coords => [$_])} ($h1)
    );
    return $mol;
  }
  
}
