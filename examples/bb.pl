use HackaMol::X::NERF;
use HackaMol;
use Modern::Perl;
use Data::Dumper;


my $bb = &init_bb;
&bb_lys($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_glu($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_thr($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_ala($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_ala($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_ala($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_lys($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_phe($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_glu($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_arg($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_gln($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_his($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_met($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_asp($bb, -120,140,180);
&extend_bb($bb,-120,140,180);
&bb_ser($bb, -120,140,180);
$bb->print_pdb;
exit;

foreach (1 .. 100){
  #&extend_bb($bb,-180,180,180);
 # &bb_his($bb);
 # &bb_lys($bb);
 # &bb_arg($bb);

 # &bb_asp($bb);
 # &bb_glu($bb);

 # &bb_ser($bb);
 # &bb_thr($bb);
 # &bb_asn($bb);
  &bb_gln($bb);

 # &bb_cys($bb);
 # &bb_secys($bb);
 # &bb_gly($bb);
 # &bb_pro($bb);

 # &bb_ala($bb);
 # &bb_val($bb);
 # &bb_ile($bb);
 # &bb_leu($bb);
 # &bb_met($bb);
 # &bb_phe($bb);
 # &bb_tyr($bb);
 # &bb_trp($bb);
  &extend_bb($bb,-60,180,180);
}
#  &extend_bb($bb, -120,140,180);
&nterm_mecap($bb);
$bb->print_xyz;

sub init_bb {
  my $nerf =  HackaMol::X::NERF->new;
  my $n  = $nerf->init(0,0,0);
  my $hm_n = HackaMol::Atom->new(symbol => 'N',name=>'N', coords=>[$n]);
  my $ca = $nerf->extend_a($n,1.45);
  my $hm_ca = HackaMol::Atom->new(symbol => 'C',name=>'CA', coords=>[$ca]);
  my $c  = $nerf->extend_ab($n,$ca, 1.52, 109.5);
  my $hm_c = HackaMol::Atom->new(symbol => 'C',name=>'C', coords=>[$c]);
  my $o  = $nerf->extend_abc($n,$ca,$c, 1.23, 120.0,0);
  my $hm_o = HackaMol::Atom->new(symbol => 'O',name=>'O', coords=>[$o]);
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

  my $n  = $nerf->extend_abc($ln,$lca,$lc, 1.33, 120.0,$psi);
  $n  = $nerf->extend_abc($ln,$lca,$lc, 1.33, 120.0,-60) if ($pro eq 'P');
  my $hm_n = HackaMol::Atom->new(symbol => 'N', name=>'N', coords=>[$n]);
  my $ca = $nerf->extend_abc($lca,$lc,$n,  1.45, 120.0,$omega);
  my $hm_ca = HackaMol::Atom->new(symbol => 'C',name=>'CA', coords=>[$ca]);
  my $c  = $nerf->extend_abc($lc,$n,$ca,  1.52, 109.5,$phi);
  my $hm_c = HackaMol::Atom->new(symbol => 'C', name=>'C',coords=>[$c]);
  my $o  = $nerf->extend_abc($n,$ca,$c,  1.23, 120,0);
  my $hm_o = HackaMol::Atom->new(symbol => 'O', name=>'O',coords=>[$o]);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $nd1    = $nerf->extend_abc($lca,$cb,$cg, 1.35, 120,-90);
  my $hm_nd1 = HackaMol::Atom->new(symbol => 'N',name=>'ND1', coords=>[$nd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$nd1, 1.34, 108.47,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C',name=>'CE1', coords=>[$ce1]);

  my $ne2    = $nerf->extend_abc($cg,$nd1,$ce1, 1.26, 107.13,0);
  my $hm_ne2 = HackaMol::Atom->new(symbol => 'N',name=>'NE2', coords=>[$ne2]);

  my $cd2    = $nerf->extend_abc($nd1,$ce1,$ne2, 1.37, 113.2,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C',name=>'CD2', coords=>[$cd2]);

  &proper_aa('HIS', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_nd1,$hm_cd2,$hm_ce1,$hm_ne2);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  &proper_aa('ALA', @atoms[-4,-3,-2,-1],$hm_cb);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg2   = $nerf->extend_abc($lc,$lca,$cb, 1.56, 109,-60);
  my $hm_cg2 = HackaMol::Atom->new(symbol => 'C',name=>'CG2', coords=>[$cg2]);

  my $og1    = $nerf->extend_abc($cg2,$lca,$cb, 1.4, 109,120);
  my $hm_og1 = HackaMol::Atom->new(symbol => 'O',name=>'OG1', coords=>[$og1]);

  &proper_aa('THR', @atoms[-4,-3,-2,-1],$hm_cb,$hm_og1,$hm_cg2);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $og    = $nerf->extend_abc($lc,$lca,$cb, 1.4, 109,120);
  my $hm_og = HackaMol::Atom->new(symbol => 'O',name=>'OG', coords=>[$og]);

  &proper_aa('SER', @atoms[-4,-3,-2,-1],$hm_cb,$hm_og);
  $bb->push_atoms($hm_cb,$hm_og);

}

sub bb_cys {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $sg    = $nerf->extend_abc($lc,$lca,$cb, 1.81, 115, 120);
  my $hm_sg = HackaMol::Atom->new(symbol => 'S',name=>'SG', coords=>[$sg]);

  &proper_aa('CYS', @atoms[-4,-3,-2,-1],$hm_cb,$hm_sg);
  $bb->push_atoms($hm_cb,$hm_sg);

}

sub bb_secys {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz;
  my $lca = $atoms[-3]->xyz;
  my $lc  = $atoms[-2]->xyz;
  my $lo  = $atoms[-1]->xyz;


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $seg    = $nerf->extend_abc($lc,$lca,$cb, 1.98, 115, 120);
  my $hm_seg = HackaMol::Atom->new(symbol => 'Se',name=>'SEG', coords=>[$seg]);

  &proper_aa('SEC', @atoms[-4,-3,-2,-1],$hm_cb,$hm_seg);
  $bb->push_atoms($hm_cb,$hm_seg);

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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG',  coords=>[$cg]);

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.41, 120,-90);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C',name=>'CD1',  coords=>[$cd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$cd1, 1.41, 120,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C',name=>'CE1',  coords=>[$ce1]);

  my $cz     = $nerf->extend_abc($cg,$cd1,$ce1, 1.41, 120,0);
  my $hm_cz  = HackaMol::Atom->new(symbol => 'C',name=>'CZ',  coords=>[$cz]);

  my $ce2    = $nerf->extend_abc($cd1,$ce1,$cz, 1.41, 120,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'C',name=>'CE2',  coords=>[$ce2]);

  my $cd2    = $nerf->extend_abc($ce1,$cz,$ce2, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C',name=>'CD2',  coords=>[$cd2]);

  &proper_aa('PHE', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd1,$hm_ce1,$hm_cz ,$hm_ce2, $hm_cd2);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.41, 120,-90);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C',name=>'CD1', coords=>[$cd1]);

  my $ce1    = $nerf->extend_abc($cb,$cg,$cd1, 1.41, 120,180);
  my $hm_ce1 = HackaMol::Atom->new(symbol => 'C',name=>'CE1', coords=>[$ce1]);

  my $cz     = $nerf->extend_abc($cg,$cd1,$ce1, 1.41, 120,0);
  my $hm_cz  = HackaMol::Atom->new(symbol => 'C',name=>'CZ', coords=>[$cz]);

  my $ce2    = $nerf->extend_abc($cd1,$ce1,$cz, 1.41, 120,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'C',name=>'CE2', coords=>[$ce2]);

  my $cd2    = $nerf->extend_abc($ce1,$cz,$ce2, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C',name=>'CD2', coords=>[$cd2]);
  
  my $oh     = $nerf->extend_abc($cd1,$ce1,$cz , 1.38, 120,180);
  my $hm_oh = HackaMol::Atom->new(symbol => 'O',name=>'OH', coords=>[$oh]);
  
  &proper_aa('TYR', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd1,$hm_ce1,$hm_cz ,$hm_ce2, $hm_cd2, $hm_oh);
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
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.36, 120,-60);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C',name=>'CD1', coords=>[$cd1]);

  my $ne1    = $nerf->extend_abc($cb,$cg,$cd1, 1.37, 108.47,180);
  my $hm_ne1 = HackaMol::Atom->new(symbol => 'N',name=>'NE1', coords=>[$ne1]);

  my $ce2    = $nerf->extend_abc($cg,$cd1,$ne1, 1.34, 111.13,0);
  my $hm_ce2 = HackaMol::Atom->new(symbol => 'N',name=>'CE2', coords=>[$ce2]);


  my $cz2    = $nerf->extend_abc($cd1,$ne1,$ce2, 1.41, 131,180);
  my $hm_cz2 = HackaMol::Atom->new(symbol => 'C',name=>'CZ2', coords=>[$cz2]);

  my $ch2    = $nerf->extend_abc($ne1,$ce2,$cz2, 1.41, 120,180);
  my $hm_ch2 = HackaMol::Atom->new(symbol => 'C',name=>'CH2', coords=>[$ch2]);

  my $cz3    = $nerf->extend_abc($ce2,$cz2,$ch2, 1.41, 120,0);
  my $hm_cz3 = HackaMol::Atom->new(symbol => 'C',name=>'CZ3', coords=>[$cz3]);

  my $ce3     = $nerf->extend_abc($cz2,$ch2,$cz3, 1.41, 120,0);
  my $hm_ce3  = HackaMol::Atom->new(symbol => 'C',name=>'CE3', coords=>[$ce3]);

  my $cd2    = $nerf->extend_abc($ch2,$cz3,$ce3, 1.41, 120,0);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C',name=>'CD2', coords=>[$cd2]);

  &proper_aa('TRP', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd1,$hm_ne1,$hm_ce2,
                  $hm_cz2,$hm_cz3,$hm_ch2,$hm_ce3,$hm_cd2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd1,$hm_ne1,$hm_ce2,
                  $hm_cz2,$hm_cz3,$hm_ch2,$hm_ce3,$hm_cd2);

}

sub bb_asp {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $od1    = $nerf->extend_abc($lca,$cb,$cg, 1.26, 120,90);
  my $hm_od1 = HackaMol::Atom->new(symbol => 'O',name=>'OD1', coords=>[$od1]);

  my $od2    = $nerf->extend_abc($od1,$cb,$cg, 1.26, 120,180);
  my $hm_od2 = HackaMol::Atom->new(symbol => 'O',name=>'OD2', coords=>[$od2]);


  &proper_aa('ASP', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_od1,$hm_od2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_od1,$hm_od2);

}


sub bb_glu {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG',coords=>[$cg]);

  my $cd   = $nerf->extend_abc($lca,$cb,$cg, 1.52, 109,180);
  my $hm_cd = HackaMol::Atom->new(symbol => 'C', name=>'CD',coords=>[$cd]);

  my $oe1    = $nerf->extend_abc($cb,$cg,$cd, 1.26, 120,90);
  my $hm_oe1 = HackaMol::Atom->new(symbol => 'O',name=>'OE1', coords=>[$oe1]);

  my $oe2    = $nerf->extend_abc($oe1,$cg,$cd, 1.26, 120,180);
  my $hm_oe2 = HackaMol::Atom->new(symbol => 'O',name=>'OE2', coords=>[$oe2]);


  &proper_aa('GLU', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd,$hm_oe1,$hm_oe2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd,$hm_oe1,$hm_oe2);

}

sub bb_lys {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $cd   = $nerf->extend_abc($lca,$cb,$cg, 1.52, 109,180);
  my $hm_cd = HackaMol::Atom->new(symbol => 'C',name=>'CD', coords=>[$cd]);

  my $ce   = $nerf->extend_abc($cb,$cg,$cd, 1.52, 109,180);
  my $hm_ce = HackaMol::Atom->new(symbol => 'C',name=>'CE', coords=>[$ce]);

  my $nz    = $nerf->extend_abc($cg,$cd,$ce, 1.5, 109,180);
  my $hm_nz = HackaMol::Atom->new(symbol => 'N',name=>'NZ', coords=>[$nz]);


  &proper_aa('LYS', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd,$hm_ce,$hm_nz);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd,$hm_ce,$hm_nz);

}

sub bb_gly {
  return;
}

sub bb_arg {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name => 'CB', coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG',coords=>[$cg]);

  my $cd   = $nerf->extend_abc($lca,$cb,$cg, 1.52, 109,180);
  my $hm_cd = HackaMol::Atom->new(symbol => 'C', name=> 'CD', coords=>[$cd]);

  my $ne    = $nerf->extend_abc($cb,$cg,$cd, 1.45, 109,180);
  my $hm_ne = HackaMol::Atom->new(symbol => 'N', name=>'NE',coords=>[$ne]);

  my $cz    = $nerf->extend_abc($cg,$cd,$ne, 1.34, 120,180);
  my $hm_cz = HackaMol::Atom->new(symbol => 'C', name=>'CZ', coords=>[$cz]);

  my $nh1    = $nerf->extend_abc($cd,$ne,$cz, 1.34, 120,180);
  my $hm_nh1 = HackaMol::Atom->new(symbol => 'n', name=>'NH1',coords=>[$nh1]);

  my $nh2    = $nerf->extend_abc($nh1,$ne,$cz, 1.34, 120,180);
  my $hm_nh2 = HackaMol::Atom->new(symbol => 'n', name=>'NH2',coords=>[$nh2]);

  &proper_aa('ARG', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd,$hm_ne,$hm_cz,$hm_nh1,$hm_nh2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd,$hm_ne,$hm_cz,$hm_nh1,$hm_nh2);

}

sub bb_val {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB', coords=>[$cb]);

  my $cg1    = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,60);
  my $hm_cg1 = HackaMol::Atom->new(symbol => 'C', name => 'CG1', coords=>[$cg1]);

  my $cg2    = $nerf->extend_abc($cg1,$lca,$cb, 1.52, 109,120);
  my $hm_cg2 = HackaMol::Atom->new(symbol => 'C', name => 'CG2', coords=>[$cg2]);

  &proper_aa('VAL', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg1,$hm_cg2);
  $bb->push_atoms($hm_cb,$hm_cg1,$hm_cg2);

}

sub bb_ile {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB', coords=>[$cb]);

  my $cg1    = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,60);
  my $hm_cg1 = HackaMol::Atom->new(symbol => 'C', name=>'CG1', coords=>[$cg1]);

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg1, 1.52, 109, 180);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C', name=>'CD1', coords=>[$cd1]);

  my $cg2    = $nerf->extend_abc($cg1,$lca,$cb, 1.52, 109,120);
  my $hm_cg2 = HackaMol::Atom->new(symbol => 'C', name=>'CG2', coords=>[$cg2]);

  &proper_aa('ILE', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg1,$hm_cg2,$hm_cd1);
  $bb->push_atoms($hm_cb,$hm_cg1,$hm_cg2,$hm_cd1);

}

sub bb_leu {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB', coords=>[$cb]);

  my $cg    = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG', coords=>[$cg]);

  my $cd1    = $nerf->extend_abc($lca,$cb,$cg, 1.52, 109, 120);
  my $hm_cd1 = HackaMol::Atom->new(symbol => 'C',name=>'CD1', coords=>[$cd1]);

  my $cd2    = $nerf->extend_abc($cd1,$cb,$cg, 1.52, 109,120);
  my $hm_cd2 = HackaMol::Atom->new(symbol => 'C',name=>'CD2', coords=>[$cd2]);

  &proper_aa('LEU', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd1,$hm_cd2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd1,$hm_cd2);

}

sub bb_asn {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB',coords=>[$cb]);

  my $cg    = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG',coords=>[$cg]);

  my $od1    = $nerf->extend_abc($lca,$cb,$cg, 1.27, 120, 90);
  my $hm_od1 = HackaMol::Atom->new(symbol => 'O', name=>'OD1',coords=>[$od1]);

  my $nd2    = $nerf->extend_abc($od1,$cb,$cg, 1.31, 120,180);
  my $hm_nd2 = HackaMol::Atom->new(symbol => 'N', name=>'ND2',coords=>[$nd2]);

  &proper_aa('ASN', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_od1,$hm_nd2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_od1,$hm_nd2);

}

sub bb_gln {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb     = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb  = HackaMol::Atom->new(symbol => 'C',name=>'CB', coords=>[$cb]);

  my $cg     = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg  = HackaMol::Atom->new(symbol => 'C',name=>'CG', coords=>[$cg]);

  my $cd    = $nerf->extend_abc($lca,$cb,$cg, 1.52, 109,180);
  my $hm_cd = HackaMol::Atom->new(symbol => 'C',name=>'CD', coords=>[$cd]);

  my $oe1    = $nerf->extend_abc($lca,$cg,$cd, 1.27, 120, 90);
  my $hm_oe1 = HackaMol::Atom->new(symbol => 'O',name=>'OE1', coords=>[$oe1]);

  my $ne2    = $nerf->extend_abc($oe1,$cg,$cd, 1.31, 120,180);
  my $hm_ne2 = HackaMol::Atom->new(symbol => 'N',name=>'NE2', coords=>[$ne2]);

  &proper_aa('GLN', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd,$hm_oe1,$hm_ne2);
  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd,$hm_oe1,$hm_ne2);


}

sub bb_pro {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($lc, $ln,$lca, 1.53, 107.,-116);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB',coords=>[$cb]);

  my $cg   = $nerf->extend_abc($ln,$lca,$cb, 1.48, 105.7,0);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG',coords=>[$cg]);

  my $cd   = $nerf->extend_abc($lca,$cb,$cg, 1.49, 110.4,0);
  my $hm_cd = HackaMol::Atom->new(symbol => 'C', name=>'CD',coords=>[$cd]);

  &proper_aa('PRO', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_cd);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_cd);

}

sub bb_met {
  my $nerf =  HackaMol::X::NERF->new;
  my $bb = shift;

  my @atoms = $bb->all_atoms;
  my $ln  = $atoms[-4]->xyz; # 1
  my $lca = $atoms[-3]->xyz; # 2
  my $lc  = $atoms[-2]->xyz; # 3 
  my $lo  = $atoms[-1]->xyz; # 4


  my $cb   = $nerf->extend_abc($ln, $lc,$lca, 1.52, 109,120);
  my $hm_cb = HackaMol::Atom->new(symbol => 'C', name=>'CB',coords=>[$cb]);

  my $cg   = $nerf->extend_abc($lc,$lca,$cb, 1.52, 109,120);
  my $hm_cg = HackaMol::Atom->new(symbol => 'C', name=>'CG',coords=>[$cg]);

  my $sd   = $nerf->extend_abc($lca,$cb,$cg, 1.81, 109,120);
  my $hm_sd = HackaMol::Atom->new(symbol => 'S', name=>'SD',coords=>[$sd]);

  my $ce   = $nerf->extend_abc($cb,$cg,$sd, 1.81, 102,90);
  my $hm_ce = HackaMol::Atom->new(symbol => 'C',name=>'CE', coords=>[$ce]);

  &proper_aa('MET', @atoms[-4,-3,-2,-1],$hm_cb,$hm_cg,$hm_sd,$hm_ce);

  $bb->push_atoms($hm_cb,$hm_cg,$hm_sd,$hm_ce);

}

sub nterm_mecap {
  my $nerf =  HackaMol::X::NERF->new;

  my $ln  = $bb->get_atoms(0)->xyz;
  my $lca = $bb->get_atoms(1)->xyz;
  my $lc  = $bb->get_atoms(2)->xyz;
  
  

  my $cme   = $nerf->extend_abc($lc,$lca,$ln, 1.5, 109,180);
  my $hm_ce = HackaMol::Atom->new(symbol => 'C', coords=>[$cme]);
#  printf ("%10.6f %10.6f %10.6f\n", @$cme); exit;
  $bb->unshift_atoms($hm_ce);
  
}

sub proper_aa{
  my $resname = shift;
  foreach my $atom (@_){
    $atom->record_name("ATOM");
    $atom->resname($resname);
  }
}



