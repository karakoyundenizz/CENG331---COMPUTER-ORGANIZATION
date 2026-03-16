#!/usr/bin/perl

# Check length of to_morse_digit_string function in .yo file
# Assumes that function starts with label "to_morse_digit_string:"
# and finishes with label "End:"

$startpos = -1;
$endpos = -1;

while (<>) {
  $line = $_;
  if ($line =~ /(0x[0-9a-fA-F]+):.* to_morse_digit_string:/) {
    $startpos = hex($1);
  }
  if ($line =~ /(0x[0-9a-fA-F]+):.* End:/) {
    $endpos = hex($1);
  }
}

if ($startpos >= 0 && $endpos > $startpos) {
  $len = $endpos - $startpos;
  print "to_morse_digit_string length = $len bytes\n";
} else {
  print "Couldn't determine to_morse_digit_string length\n";
}
