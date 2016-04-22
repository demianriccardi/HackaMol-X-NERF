use HackaMol::X::NERF;
use HackaMol;
use Modern::Perl;
use Data::Dumper;

my $bb = &init_bb;
foreach (1 .. 100){
  &extend_bb($bb);
  &bb_his($bb);
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
  my $pro = shift || 'np';
  
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
  $c  = $nerf->extend_abc($lc,$n,$ca,  1.52, 109.5,-60) if ($pro eq 'P');
  my $hm_c = HackaMol::Atom->new(symbol => 'C', coords=>[$c]);
  my $o  = $nerf->extend_abc($n,$ca,$c,  1.23, 120,0);
  my $hm_o = HackaMol::Atom->new(symbol => 'O', coords=>[$o]);
  my $residue = HackaMol::AtomGroup->new(atoms=>[$hm_n,$hm_ca,$hm_c,$hm_o]);
  $bb->push_groups($residue); 
  
}

sub bb_his {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', coords=>[$cg]);

  my $nd1    = $nerf->extend_abc($lca,$cb,$cg, 1.35, 120,-90);
  my $hm_nd1 = HackaMol::Atom->new(symbol => 'N', coords=>[$nd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$nd1, 1.34, 108.47,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'N', coords=>[$ce1]);

  my $ne2    = $nerf->extend_abc($cg,$nd1,$ce1, 1.26, 107.13,0);
  my $hm_ne2 = HackaMol::Atom->new(symbol => 'N', coords=>[$ne2]);

  my $cd2    = $nerf->extend_abc($nd1,$ce1,$ne2, 1.37, 113.2,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd2]);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_nd1,$hm_cd2,$hm_ce1,$hm_ne2); 

}
