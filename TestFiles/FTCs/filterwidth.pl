#! /usr/bin/perl

use PDL;
use PGPLOT;
use PDL::Graphics::PGPLOT;
use PDL::AutoLoader;
use PDL::Math;

dev '/xs';
#dev '/vcps';
plotting_setup (5,1);


#$filterfile = "f814w_sys.ftc";
$filterfile = "J.ftc";

($w, $t) = rcols $filterfile; 

line $w, $t; hold; 
$hh = 0.5*max $t;  # <- the height at half maximum transmission;

$aha = abs($t-$hh)-($t-$hh);

($whh, $aha) = where ($w,$aha, $aha == 0);

line $whh, $aha, {color=>red};

($minw, $maxw) = minmax $whh; 

line (pdl ($minw, $maxw), $hh*ones(2),{color=>red});
print "$minw $maxw \n";
 

