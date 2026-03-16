















void touch2(unsigned int val1)
{
    vlevel = 2; // Part of the validation protocol
    if (val1 == cookie) {
    	printf("Touch2!: You called touch2(0x%.8x)\n", val1);
    	validate(2);
    } else {
    	printf("Misfire: You called touch2(0x%.8x)\n", val1);
    	fail(2);
    }
    exit(0);
}


















