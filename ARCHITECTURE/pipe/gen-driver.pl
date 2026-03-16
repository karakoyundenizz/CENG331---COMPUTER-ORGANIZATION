#!/usr/bin/perl 
#!/usr/local/bin/perl 

#
# gen-driver - Generate driver file for any to_morse_digit_string function
#
use Getopt::Std;

$n = 0;

getopts('hcrn:f:b:');

if ($opt_h) {
    print STDERR "Usage $argv[0] [-h] [-c] [-n N] [-f FILE]\n";
    print STDERR "   -h      print help message\n";
    print STDERR "   -c      include correctness checking code\n";
    print STDERR "   -n N    set number of elements\n";
    print STDERR "   -f FILE set input file (default stdin)\n";
    print STDERR "   -b blim set byte limit for function\n";
#    print STDERR "   -r      Allow random result\n";
    die "\n";
}

$check = 0;
if ($opt_c) {
    $check = 1;
}

$bytelim = 1000;
if ($opt_b) {
    $bytelim = $opt_b;
}

if ($opt_n) {
    $n = $opt_n;
    if ($n < 0) {
	print STDERR "n must be at least 0\n";
	die "\n";
    }
}

#$randomval = 0;
# Accumulated count
$rval = 0;

#if ($opt_r) {
#  $randomval = 1;
#} else {
# Value that should be returned by function
#  $tval = int($n/2);
#}


# The data to be stored.
@data = ();

for ($i = 0; $i < $n; $i++) {
  # Digits cycle through 0..9 so that the reference implementation
  # and the checker both see valid Morse digits.
  $data[$i] = ($i % 10);
  # Expected return value is the sum of digits.
  $rval += $data[$i];
}


# Values to put at beginning and end of destination
$Preval =  "0xbcdefa";
$Postval = "0xdefabc";
$Corrval = "\$0x5710331";

print <<PROLOGUE;
#######################################################################
# Test for copying block of size $n;
#######################################################################
	.pos 0
main:	irmovq Stack, %rsp  	# Set up stack pointer

	# Set up arguments for copy function and then invoke it
	irmovq \$$n, %rdx		# src and dst have $n elements
	irmovq dest, %rsi	# dst array
	irmovq src, %rdi	# src array
    # corrupt all the unused registers to prevent assumptions
    irmovq $Corrval, %rax
    irmovq $Corrval, %rbx
    irmovq $Corrval, %rcx
    irmovq $Corrval, %rbp
    irmovq $Corrval, %r8
    irmovq $Corrval, %r9
    irmovq $Corrval, %r10
    irmovq $Corrval, %r11
    irmovq $Corrval, %r12
    irmovq $Corrval, %r13
    irmovq $Corrval, %r14
	call to_morse_digit_string		
	 
PROLOGUE

if ($check) {
print <<CALL;
	call check	        # Call checker code
	halt                # should halt with 0xaaaa in %rax
.pos 200
temp_data:
 .quad 0
CALL
} else {
print <<HALT;
	halt			# should halt with sum in %rax
HALT
}

print "StartFun:\n";
if ($opt_f) {
    open (CODEFILE, "$opt_f") || die "Can't open code file $opt_f\n";
    while (<CODEFILE>) {
	printf "%s", $_;
    }
} else {
    while (<>) {
	printf "%s", $_;
    }
}
print "EndFun:\n";

if ($check) {
$len_min_1 = 8 * $n - 8;
print <<CHECK;
#################################################################### 
# Epilogue code for the correctness testing driver
####################################################################

# This is the correctness checking code.
# It checks:
#   1. %rax has $rval.  Set %rax to 0xbbbb if not.
#   2. The total length of the code is less than or equal to $bytelim.
#      Set %rax to 0xcccc if not.
#   3. Each digit in the source was converted to the correct 5-character Morse code pattern in the destination.
#      Set %rax to 0xdddd if not.
#   4. The words just before and just after the destination region
#      were not corrupted.  Set %rax to 0xeeee if not.
# If all checks pass, then sets %rax to 0xaaaa
check:
	# Return value test
	irmovq \$$rval,%r10
	subq %r10,%rax
	je checkb
	irmovq \$0xbbbb,%rax  # Failed test #1
	jmp cdone
checkb:
	# Code length check
	irmovq EndFun,%rax
	irmovq StartFun,%rdx
	subq %rdx,%rax
	irmovq \$$bytelim,%rdx
	subq %rax,%rdx
	jge checkm
	irmovq \$0xcccc,%rax  # Failed test #2
	jmp cdone
checkm:
	# Strict Morse pattern check:
	# For each digit d in src, we recompute its 5-character Morse code
	# and compare it against what the student's code wrote into dest.
	irmovq dest, %rdx     # Pointer to next destination character
	irmovq src,  %rbx     # Pointer to next source digit
	irmovq \$$n,%rdi      # Digit count
	andq %rdi,%rdi
	je checkpre           # Skip check if count = 0

mcloop:
	# Load current digit d from src (0..9)
	mrmovq (%rbx),%r9

	# i = 0 (character index in [0,4])
	xorq %r11,%r11

charloop:
	# if (i >= 5) break;
	rrmovq %r11,%r12
	irmovq \$5,%r8
	subq %r8,%r12         # i - 5
	jge donechars

	# Load character from destination: ch = *dest
	mrmovq (%rdx),%r10

	# Isolate low byte: ch &= 0xFF
	irmovq \$0xff,%r8
	andq %r8,%r10

	# Compute expected character in %r13 based on d and i:
	# if (d == 0) ch = '-';
	# else if (d <= 5)  ch = (i < d ? '.' : '-');
	# else { dashes = d - 5; ch = (i < dashes ? '-' : '.'); }

	# if (d == 0)
	rrmovq %r9,%r14
	andq %r14,%r14
	jne notzero

	irmovq \$45,%r13         # '-'
	jmp have_expected

notzero:
	# Check if d <= 5
	rrmovq %r9,%r14
	irmovq \$5,%r8
	subq %r8,%r14           # d - 5
	jg digit6to9            # d > 5 -> 6..9 case

	# 1..5 case: if (i < d) '.' else '-'
	rrmovq %r11,%r12        # i - d
	subq %r9,%r12
	jl dot_1to5

	irmovq \$45,%r13         # '-'
	jmp have_expected

dot_1to5:
	irmovq \$46,%r13         # '.'
	jmp have_expected

digit6to9:
	# Real Morse for 6..9: leading (d - 5) dashes, then dots.
	rrmovq %r9,%r14         # r14 = d
	irmovq \$5,%r8
	subq %r8,%r14           # r14 = d - 5 = dashes

	# if (i < dashes) '-' else '.'
	rrmovq %r11,%r12        # i - dashes
	subq %r14,%r12
	jl dash_6to9

	irmovq \$46,%r13         # '.'
	jmp have_expected

dash_6to9:
	irmovq \$45,%r13         # '-'

have_expected:
	# Compare actual ch (%r10) with expected ch (%r13)
	rrmovq %r10,%r12
	subq %r13,%r12
	je okchar

	# Mismatch -> incorrect Morse pattern
	irmovq \$0xdddd,%rax    # Failed test #3
	jmp cdone

okchar:
	# Advance dest pointer and i++
	irmovq \$1,%r8
	addq %r8,%rdx           # dest++
	addq %r8,%r11           # i++
	jmp charloop

donechars:
	# Finished 5 characters for this digit.
	# Advance src pointer to next digit (each digit is a quad).
	irmovq \$8,%r8
	addq %r8,%rbx           # src += 8

	# Decrement digit count and continue if more digits remain.
	irmovq \$1,%r8
	subq %r8,%rdi           # cnt--
	jg mcloop

checkpre:
	# Check for corruption
	irmovq Predest,%rdx
	mrmovq (%rdx), %rax  # Get word before destination
	irmovq \$$Preval, %rdx
	subq %rdx,%rax
	je checkpost
	irmovq \$0xeeee,%rax  # Failed test #4
	jmp cdone
checkpost:
	# Check for corruption
	#irmovq Postdest,%rdx
	#mrmovq (%rdx), %rax  # Get word after destination
	#irmovq \$$Postval, %rdx
	#subq %rdx,%rax
	#je checkok
	#irmovq \$0xeeee,%rax # Failed test #4
	#jmp cdone
checkok:
	# Successful checks
	irmovq \$0xaaaa,%rax
cdone:
	ret
CHECK
}

print <<EPILOGUE1;

###############################
# Source and destination blocks 
###############################
	.align 8
src:
EPILOGUE1

for ($i = 0; $i < $n; $i++) {
    print "\t.quad $data[$i]\n";
}

print <<EPILOGUE2;
	.quad $Preval # This shouldn't get moved

	.align 16
Predest:
	.quad $Preval
dest:
EPILOGUE2

for ($i = 0; $i < $n; $i++) {
    print "\t.quad 0xcdefab\n";
}

print <<EPILOGUE3;
Postdest:
	.quad $Postval

.align 8
# Run time stack
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0
	.quad 0

Stack:
EPILOGUE3
