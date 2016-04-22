use HackaMol::X::NERF;
use HackaMol;
use Modern::Perl;
use Data::Dumper;

my $bb = &init_bb;
foreach (1 .. 100){
  &extend_bb($bb,-180,180,180);
  #&extend_bb($bb,-120,140,180);
 # &bb_his($bb);
 # &bb_thr($bb);
 # &bb_ser($bb);
 # &bb_ala($bb);
 # &bb_phe($bb);
 # &bb_tyr($bb);
  &bb_trp($bb);
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
  my ($phi,$psi,$omega) = @_;
  ($phi,$psi,$omega) = (180,180,180) unless @_;
  my $pro = shift || 'np';
  
  my $last_group = $bb->count_groups-1;
  my @atoms = $bb->get_groups($last_group)->all_atoms;
  my $ln  = $atoms[0]->xyz;
  my $lca = $atoms[1]->xyz;
  my $lc  = $atoms[2]->xyz;
  my $lo  = $atoms[3]->xyz;

  my $n  = $nerf->extend_abc($lca,$lo,$lc, 1.33, 120.0,$omega);
  my $hm_n = HackaMol::Atom->new(symbol => 'N', coords=>[$n]);
  my $ca = $nerf->extend_abc($lca,$lc,$n,  1.45, 120.0,$phi);
  my $hm_ca = HackaMol::Atom->new(symbol => 'C', coords=>[$ca]);
  my $c  = $nerf->extend_abc($lc,$n,$ca,  1.52, 109.5,$psi);
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
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C', coords=>[$ce1]);

  my $ne2    = $nerf->extend_abc($cg,$nd1,$ce1, 1.26, 107.13,0);
  my $hm_ne2 = HackaMol::Atom->new(symbol => 'N', coords=>[$ne2]);

  my $cd2    = $nerf->extend_abc($nd1,$ce1,$ne2, 1.37, 113.2,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd2]);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_nd1,$hm_cd2,$hm_ce1,$hm_ne2); 

}

sub bb_ala {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', coords=>[$cb]);

  $bb->push_atoms($hm_cb); 

}

sub bb_thr {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', coords=>[$cb]);

  my $cg2   = $nerf->extend_abc($lc,$lca,$cb, 1.56, 109,-60);
  my $hm_cg2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cg2]);

  my $og1    = $nerf->extend_abc($cg2,$lca,$cb, 1.4, 109,120);
  my $hm_og1 = HackaMol::Atom->new(symbol => 'O', coords=>[$og1]);

  $bb->push_atoms($hm_cb,$hm_og1,$hm_cg2); 

}
    
sub bb_ser {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', coords=>[$cb]);

  my $og1    = $nerf->extend_abc($lc,$lca,$cb, 1.4, 109,120);
  my $hm_og1 = HackaMol::Atom->new(symbol => 'O', coords=>[$og1]);

  $bb->push_atoms($hm_cb,$hm_og1);

}

sub bb_phe {
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

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.41, 120,-90);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$cd1, 1.41, 120,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C', coords=>[$ce1]);

  my $cz     = $nerf->extend_abc($cg,$cd1,$ce1, 1.41, 120,0);
  my $hm_cz  = HackaMol::Atom->new(symbol => 'C', coords=>[$cz]);

  my $ce2    = $nerf->extend_abc($cd1,$ce1,$cz, 1.41, 120,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'C', coords=>[$ce2]);

  my $cd2    = $nerf->extend_abc($ce1,$cz,$ce2, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd2]);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd1,$hm_ce1,$hm_cz ,$hm_ce2, $hm_cd2);

}

sub bb_tyr {
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

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.41, 120,-90);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$cd1, 1.41, 120,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C', coords=>[$ce1]);

  my $cz     = $nerf->extend_abc($cg,$cd1,$ce1, 1.41, 120,0);
  my $hm_cz  = HackaMol::Atom->new(symbol => 'C', coords=>[$cz]);

  my $ce2    = $nerf->extend_abc($cd1,$ce1,$cz, 1.41, 120,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'C', coords=>[$ce2]);

  my $cd2    = $nerf->extend_abc($ce1,$cz,$ce2, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd2]);
  
  my $oh     = $nerf->extend_abc($cd1,$ce1,$cz , 1.38, 120,180);
  my $hm_oh = HackaMol::Atom->new(symbol => 'O', coords=>[$oh]);
  

  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd1,$hm_ce1,$hm_cz ,$hm_ce2, $hm_cd2, $hm_oh);

}

sub bb_trp {
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

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.36, 120,-60);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd1]);

  my $ne1    = $nerf->extend_abc($cb,$cg,$cd1, 1.37, 108.47,180);
  my $hm_ne1 = HackaMol::Atom->new(symbol => 'N', coords=>[$ne1]);

  my $ce2    = $nerf->extend_abc($cg,$cd1,$ne1, 1.34, 111.13,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'N', coords=>[$ce2]);


  my $cz2    = $nerf->extend_abc($cd1,$ne1,$ce2, 1.41, 131,180);
  my $hm_cz2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cz2]);

  my $ch2    = $nerf->extend_abc($ne1,$ce2,$cz2, 1.41, 120,180);
  my $hm_ch2 = HackaMol::Atom->new(symbol => 'C', coords=>[$ch2]);

  my $cz3    = $nerf->extend_abc($ce2,$cz2,$ch2, 1.41, 120,0);
  my $hm_cz3 = HackaMol::Atom->new(symbol => 'C', coords=>[$cz3]);

  my $ce3     = $nerf->extend_abc($cz2,$ch2,$cz3, 1.41, 120,0);
  my $hm_ce3  = HackaMol::Atom->new(symbol => 'C', coords=>[$ce3]);

  my $cd2    = $nerf->extend_abc($ch2,$cz3,$ce3, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C', coords=>[$cd2]);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd1,$hm_ne1,$hm_ce2,
                  $hm_cz2,$hm_cz3,$hm_ch2,$hm_ce3,$hm_cd2);

}
