#Made by Mathias Mangino
#Send on 10-05-2021

use strict;

my $f1 = shift or die "multifasta file missing\n";
my $f2 = shift or die "missing coordinate file\n";
my $f3 = shift;
my $f4 = shift;

# por defecto convertir to_uppp si recibo flag no

extractor($f1,$f2,$f3,$f4);

 
sub extractor {
	my $archivo_multi_fasta = $_[0];
	my $archivo_con_coordenadas = $_[1];
	my $out_put = $_[2];
	my $flag_not_to_upper = $_[3];
		
	my %hash_secuencias = leer_multi_fasta($archivo_multi_fasta);
	
	open (READ,$archivo_con_coordenadas) or die "no file exists $archivo_con_coordenadas\n";
	while (<READ>) {
		my ($nombre, $coor_incio, $coor_fin, $nombre_prot, @resto) =  split(/\s+/, $_);
		my $contig = $hash_secuencias{$nombre};
		my $secuencia = extraer_secuencia($contig,$coor_incio,$coor_fin);
		
		formatear_secuencia(\$secuencia);
		
		unless($flag_not_to_upper){
			$secuencia = uc($secuencia);
		}
		
		my $fasta_string = ">$nombre"."_$nombre_prot [$coor_incio $coor_fin] ".join(" ",@resto)."\n".$secuencia."\n";;
		
		if ($out_put){
			open(FILE,">>./$out_put");
			print FILE "$fasta_string";
			close(FILE);
		}else{
			print "$fasta_string";
		}
	}
	close(READ);
}

sub formatear_secuencia{
	my $secuencia = $_[0];
	$$secuencia =~ s/(.{80})/$1\n/g;
	$$secuencia =~ s/(.{10})/$1 /g;
}

sub leer_multi_fasta {
	my $archivo_multi_fasta = $_[0];
	my %hash_retorno = ();
	my $secuencias = '';
	my $nombre_secuencias = '';
	
	open (READ,$archivo_multi_fasta) or die "no file exists $archivo_multi_fasta\n";
	while(my $line = <READ>){

		chomp($line);
		
		if($line =~ />(.*)/) {
        		if($secuencias) { $hash_retorno{$nombre_secuencias} = $secuencias ; }
        		$nombre_secuencias = $1;
			$nombre_secuencias =~ s/ .*//;
        		$secuencias = '';
		}else {
			$line =~ s/\s+//g;
        		$secuencias .= $line;
   		}
	}
	$hash_retorno{$nombre_secuencias} = $secuencias;
	close(READ);
	
	return %hash_retorno;
}

sub extraer_secuencia {
	my $secuencia = $_[0];
	my $coordenada_inicio = $_[1];
	my $coordenada_fin = $_[2];
#	if ($coordenada_inicio == 0){ die "la coordenada de inicio no puede tomar valores menores a uno";}
	
	$coordenada_inicio--;
	$coordenada_fin--;
	my $largo = abs ($coordenada_fin-$coordenada_inicio)+1;
	my $retorno = '';
 	
	if ($coordenada_inicio <= $coordenada_fin) {
		$retorno = substr($secuencia,$coordenada_inicio,$largo);
	}else{
		my $seq_aux = substr($secuencia,$coordenada_fin,$largo);
		$retorno = inversa_complementaria(\$seq_aux);
	}
	return $retorno;
}

sub inversa_complementaria {
	my $dna_seq_ref = $_[0];
	my $dna_seq_rev_comp = reverse $$dna_seq_ref;
	$dna_seq_rev_comp =~ tr/ACGTacgt/TGCAtgca/;
	return $dna_seq_rev_comp;
}
