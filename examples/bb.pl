use HackaMol::X::NERF;
use HackaMol;
use Modern::Perl;
use Data::Dumper;

my $bb = &init_bb;
foreach (1 .. 100){
  &extend_bb($bb);
}
$bb->print_xyz;

sub init_bb {
  my $nerf =  HackaMol::X::NERF->new;
  my $n  = $nerf->init(0,0,0);
  my $hm_n = HackaMol::Atom->new(symbol => 'N', coords=>[$n]);
  my $ca = $nerf->extend_a($n,1.45);
  my $hm_ca = HackaMol::Atom->new(symbol => 'C', coords=>[$ca]);
  my $c  = $nerf->extend_ab($n,$ca, 1.52, 109.5);
  my $hm_c = HackaMol::Atom->new(symbol => 'C', coords=>[$c]);
  my $o  = $nerf->extend_abc($n,$ca,$c, 1.23, 120.0,0);
  my $hm_o = HackaMol::Atom->new(symbol => 'O', coords=>[$o]);
  my $residue = HackaMol::AtomGroup->new(atoms=>[$hm_n,$hm_ca,$hm_c,$hm_o]);
  return (HackaMol::Molecule->new(groups=>[$residue]));
}

sub extend_bb {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;
  my $pro = shift;
  
  my $last_group = $bb->count_groups-1;
  my @atoms = $bb->get_groups($last_group)->all_atoms;
  my $ln  = $atoms[0]->xyz;
  my $lca = $atoms[1]->xyz;
  my $lc  = $atoms[2]->xyz;
  my $lo  = $atoms[3]->xyz;

  my $n  = $nerf->extend_abc($lca,$lo,$lc, 1.33, 120.0,180);
  my $hm_n = HackaMol::Atom->new(symbol => 'N', coords=>[$n]);
  my $ca = $nerf->extend_abc($lca,$lc,$n,  1.45, 120.0,180);
  my $hm_ca = HackaMol::Atom->new(symbol => 'C', coords=>[$ca]);
  my $c  = $nerf->extend_abc($lc,$n,$ca,  1.52, 109.5,180);
  $c  = $nerf->extend_abc($lc,$n,$ca,  1.52, 109.5,-60) if ($pro);
  my $hm_c = HackaMol::Atom->new(symbol => 'C', coords=>[$c]);
  my $o  = $nerf->extend_abc($n,$ca,$c,  1.23, 120,0);
  my $hm_o = HackaMol::Atom->new(symbol => 'O', coords=>[$o]);
  my $residue = HackaMol::AtomGroup->new(atoms=>[$hm_n,$hm_ca,$hm_c,$hm_o]);
  $bb->push_groups($residue); 
  
}

sub bb_his {
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[0]->xyz;
  my $lca = $atoms[1]->xyz;
  my $lc  = $atoms[2]->xyz;
  my $lo  = $atoms[3]->xyz;

  my @zmat = map{ 
    sprintf("%s 0 %10.6f %10.6f %10.6f", $_->symbol,@{$_->xyz})
  } @atoms;
  my $zmat = join("\n",@zmat); 

  $zmat = $zmat . 
'
  C 2 CC 3 tet 4 75   
  C 5 CC 2 tet 1 135
  N 6 CN 5 trig  2 87
  C 7 CN 6 ring  5 180
  C 6 CdC 7 ring  2 115 
  N 9 CN 6 ring 5 175 


  CdC = 1.34
  CC = 1.54
  CN = 1.31
  trig = 120
  tet = 109
  ring = 106
';
  my $mol = HackaMol->new->read_string_mol($zmat,'zmat');
  $mol->print_xyz; exit;
  #$bb->push_atoms($mol->get_atoms($_)) foreach 4 .. $mol->natoms -1;  

}
