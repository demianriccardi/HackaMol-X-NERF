use Modern::Perl;
use Data::Dumper;
use WebService::Plotly;
 
my $plotly = WebService::Plotly->new( un  => 'demianriccardi', 
                                      key => 'ozdxs0hpwe');
 
 

my @ss = map {[split]}(
   '1.9520     0.0005',
   '1.9560     0.0007',
   '1.9600     0.0004',
   '1.9640     0.0009',
   '1.9680     0.0016',
   '1.9720     0.0007',
   '1.9760     0.0013',
   '1.9800     0.0014',
   '1.9840     0.0018',
   '1.9880     0.0026',
   '1.9920     0.0026',
   '1.9960     0.0037',
   '2.0000     0.0050',
   '2.0040     0.0097',
   '2.0080     0.0110',
   '2.0120     0.0125',
   '2.0160     0.0187',
   '2.0200     0.0320',
   '2.0240     0.0437',
   '2.0280     0.0801',
   '2.0320     0.1129',
   '2.0360     0.1029',
   '2.0400     0.0845',
   '2.0440     0.0685',
   '2.0480     0.0580',
   '2.0520     0.0454',
   '2.0560     0.0342',
   '2.0600     0.0295',
   '2.0640     0.0297',
   '2.0680     0.0200',
   '2.0720     0.0207',
   '2.0760     0.0185',
   '2.0800     0.0135',
   '2.0840     0.0129',
   '2.0880     0.0126',
   '2.0920     0.0097',
   '2.0960     0.0095',
   '2.1000     0.0079',
   '2.1040     0.0066',
   '2.1080     0.0058',
   '2.1120     0.0047',
   '2.1160     0.0045',
   '2.1200     0.0058',
   '2.1240     0.0045',
   '2.1280     0.0036',
   '2.1320     0.0029',
   '2.1360     0.0024',
   '2.1400     0.0032',
   '2.1440     0.0020',
   '2.1480     0.0014',
);


my @ds = map{$_->[0]} @ss;
my @ps = map{$_->[1]} @ss;
my $response = $plotly->plot( \@ds, \@ps );
print Dumper $response;

print "url is: $response->{url} \n";
print "filename on our server is: $response->{filename} \n";
