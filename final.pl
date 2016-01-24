# January 2016
# Julia Raices
# Program to do everything I ever did with this tables
#!/usr/bin/perl

use strict; # doesn't let you use variables without declaring them

# declared variables:
my ($exps, $dnds, $age, $lines, $line, $lin, $i, $j, $k, $l, $m, $n, $o, $key); # counters and strings
my ($Group, $XorA, $OorN); # strings to be printed
my (@exp, @ages, @dndss); # arrays with raw data
my (%exp_data, %dnds_data, %age_data); # hashs with printable data

# undefine variables, so that they don't - by chance - start not empty
undef $exps;
undef $dnds;
undef $age;
undef $lines;
undef $line;
undef $lin;
undef $i;
undef $j;
undef $k;
undef $l;
undef $m;
undef $n;
undef $o;
undef $key;
undef $Group;
undef $XorA;
undef $OorN;
undef @exp;
undef @ages;
undef @dndss;
undef %exp_data;
undef %dnds_data;
undef %age_data;

# sees what is the argument (after the program was called) and stores it. It should be the names of files to be used by the program
for($i=0; $i<=$#ARGV; $i+=2){
	if($ARGV[$i] eq "-exp"){
		$exps=$ARGV[$i+1];
	}
	if($ARGV[$i] eq "-dnds"){
		$dnds=$ARGV[$i+1];
	}
	if($ARGV[$i] eq "-age"){
		$age=$ARGV[$i+1];
	}
}

# check if files have been designated correctly, if they are not, it tells you how to use the program and exits
if($exps eq "" || $dnds eq "" || $age eq ""){
	print STDERR "Program usage: perl final.pl -exp EXPRESSION_TABLE -dnds DNDS_TABLE -age AGE_TABLE\n\t-exp: Table with genes' expression during each spermatogenesis phase.\n\t-dnds: Table with genes' data of dN, dS, pN, pS, etc.\n\t-age: Table with genes' age data.\n";
	exit 0;
}

# open log file and adds data from this time program was run
open (LOG, ">>final.log");
print LOG "\n".(localtime)."\nTable for expression data: $exps\nTable for dN/dS data: $dnds\nTable for age data: $age\nOutput file: final.output\n";

# finds if input files actually exists and are openable, if not, prints error in log and in stderr and exists program
unless(open(EXP, $exps) && open(DNDS, $dnds) && open(AGE, $age)){
	print STDERR "Couldn't open files, please check if files exist and if you have permission to read them.\n";
	print LOG "Error opening files.\n";
	exit 0;
}

# open output files, and prints header in it
open(OUTPUT, ">final.output");
print OUTPUT "id\tMitosis\tMeiosis\tPostMeiosis\tMitosisPostmeiosis\tMeiosisPostmeiosis\tMisosisMeiosis\tXorA\tGroup\tage\tbias\tdn\tds\tdnds\tpn\tps\tpnps\tln\tls\tfetValue\tneutrality\talpha\n";

# get dn and ds data from dNdS table
while(<DNDS>){
	$line=$_;
	chomp $line;
	@dndss=split(/\t/, $line);
	if($line=~/Symbol/){}
	$dnds_data{$dndss[13]}="$dndss[4]\t$dndss[2]\t$dndss[4]/$dndss[2]\t$dndss[5]\t$dndss[3]\t$dndss[5]/$dndss[3]\t$dndss[6]\t$dndss[7]\t$dndss[8]\t$dndss[10]\t$dndss[11]";
	$k++;
	#print "$dndss[13]\n";
}

# get age and bias data from age table
while(<AGE>){
	$lin=$_;
	chomp $lin;
	@ages=split(/\t/, $lin);
	if($lin=~/Symbol/){}
	if($ages[2]==0){
		$OorN="old";
	}
	else{
		$OorN="new";
	}
	$age_data{$ages[0]}="$OorN\t$ages[3]";
	$l++;
	#print "$ages[0]\n";
}

# getting id, branch, expression in all spermatogenesis phases, comparisson of expression in each 2 phases, group of each gene from expression data table
while(<EXP>){
	$lines=$_;
	chomp $lines;
	@exp=split(/\|/, $lines);
	if($lines=~/Symbol/){}
	# make new group for genes according to their classes, and adress if they are Autosomal or X-linked
	if($exp[3] eq "---"){
		$exp[3] = $exp[1];
	}
	if($exp[4] eq "arm_X"){
		$XorA="X";
	}
	else{
		$XorA="A";
	}
	if($exp[11] eq "1" || $exp[11] eq "2" || $exp[11] eq "3"){
		$Group="PostMeiotic";
	}
	elsif($exp[11] eq "4" || $exp[11] eq "7" || $exp[11] eq "12"){
		$Group="Meiotic";
	}
	elsif($exp[11] eq "5"){
		$Group="MeioticPostmeiotic";
	}
	elsif($exp[11] eq "6" || $exp[11] eq "8" || $exp[11] eq "10"){
		$Group="Mitotic";
	}
	elsif($exp[11] eq "9"){
		$Group="MitoticMeiotic";
	}
	elsif($exp[11] eq "11"){
		$Group="TheV";
	}
	elsif($exp[11] eq "13"){
		$Group="Equal";
	}
	else{
		$Group="Impossible";
	}
	$exp_data{$exp[3]}="$exp[3]\t$exp[5]\t$exp[6]\t$exp[7]\t$exp[8]\t$exp[9]\t$exp[10]\t$XorA\t$Group";
	$j++;
	#print "$exp[3]\n";
}
#$key = $exp[3];
# prints data in output. gives an NA (Not Avaiable) for not avaiable dn/ds data.
foreach $key (keys %exp_data){
#print "in foreach\n";
	#print "got above if\n";
	if($age_data{$exp[3]}){
		$m++;
		if($dnds_data{$key} ne ""){
			print OUTPUT "$exp_data{$key}\t$age_data{$key}\t$dnds_data{$key}\n";
			$n++;
		}
		else{
			print OUTPUT "$exp_data{$key}\t$age_data{$key}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
			$o++;
		}
	}
}

# print what was done and found in log file
print LOG "\tThere were $j genes in Expression Table.\n\t\t$l genes in Age Table.\n\t\t$k genes in dNdS Table.\n\t$m genes were found in both the Age and Expression Table.\n\t\t$n genes were printed with data from all three tables.\n\t\t$o genes were printed with data only for the Expression and Age tables.\n";

# close files used and exit program.
close LOG;
close OUTPUT;
close EXP;
close AGE;
close DNDS;

exit;

