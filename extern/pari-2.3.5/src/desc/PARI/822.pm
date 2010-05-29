#!/usr/bin/perl -w
#Copyright (C) 2003  The PARI group.
#
#This file is part of the GP2C package.
#
#PARI/GP is free software; you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software
#Foundation. It is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY WHATSOEVER.
#
#Check the License for details. You should have received a copy of it, along
#with the package; see the file 'COPYING'. If not, write to the Free Software
#Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

#Based on Debconf::Format::822 by Joey Hess <joey@kitenet.net>.

package PARI::822;
use strict;

sub new { bless {} }

=d1 NAME

PARI::822::read -- Read Description files.

=head1 SYNOPSIS

$database->PARI::822::new();
$database->read($filename,$mode)

PARI::822::read(\%database,$filename,$mode)

=head1 DESCRIPTION

read the database file $filename and merge the information in the database.

Mode is a bitmap flag
mode&1: new values cannot overwrite old ones.
mode&2: new functions are not allowed.

=cut

sub read 
{
	
	local $/="\n";
        local *FILE;
	
	my ($ret,$file,$check)=@_;
        $check=0 if (!defined($check));
	my $invars=0;
	my ($key, $value);
	my $line;
        my $entry;

        open FILE,"<$file";
	while ($line = <FILE>) 
        {
		chomp $line;
		if ($invars) 
                {
			if ($line =~ /^\s/) 
                        {
				$line =~ s/^\s+//;
				$value.="\n$line";
				next;
			}
			else
                        {
				$invars=0;
                                $value =~ s/\s*$//;
                                if (!defined($ret->{$entry}->{$key}))
                                {
                                        $ret->{$entry}->{$key}=$value;
                                }
                                elsif (($check&1) and $ret->{$entry}->{$key} ne $value)
                                {
                                        die ("Unmatched data: $entry: $key: $ret->{$entry}->{$key} ne $value");
                                }
			}
		}
                next if ($line eq '');
		($key, $value)=split(/:\s*/, $line, 2);
                die("Bad entry $key") if (!defined($value));
		if ($key eq 'Function') 
                {
                        $entry=$value;
		}
                die("New function $value") if (($check&2) and !defined($ret->{$entry}));
                $invars=1;
	}
        die("$file is not terminated") if ($invars);
        return 0;
}

=d1 NAME

PARI::822::write -- Write Description files.

=head1 SYNOPSIS

$database->PARI::822::new();
$database->write($filename)

PARI::822::write(\%database,STREAM)


=head1 DESCRIPTION

output a database to STREAM in canonical 822 format.

=cut


sub write 
{
        my @order=("Function","Class","Section","C-Name","Prototype","Type","Help","Description");
        my %knowfields=map {$_ => 1}  @order;
	my %data=%{shift()};
        my $STREAM=shift;
        defined($STREAM) or $STREAM=*STDOUT;
	foreach my $func (sort keys %data) 
        {
	        foreach my $field (@order)
                {
		        my $val=$data{$func}->{$field};
                        next if (!defined($val));
                        $val =~ s/\n/\n /g;
		        print $STREAM $field.": $val\n";
                }
	        foreach my $field (sort keys %{$data{$func}})
                {
                        next if ($knowfields{$field});
		        my $val=$data{$func}->{$field};
                        $val =~ s/\n/\n /g;
		        print $STREAM $field.": $val\n";
                }
                print $STREAM "\n";
	}
}
1;
