typedef long word_t;
/*
 * arr: length len, each arr[i] is a digit 0..9
 * buff: for each digit we write 5 chars (no '\0')
 * return: sum of digits
 *
 * It produces the following patterns:
 *   0 → "-----"
 *   1 → ".----"
 *   2 → "..---"
 *   3 → "...--"
 *   4 → "....-"
 *   5 → "....."
 *   6 → "-...."
 *   7 → "--..."
 *   8 → "---.."
 *   9 → "----."
 */
word_t to_morse_digit_string(word_t *arr, unsigned char *buff, word_t len)
{
    word_t sum = 0;

    while (len > 0) {
        word_t d = *arr++;
        sum += d;

        unsigned char *temp = buff;

        if (d == 0) {
            // 0 → "-----"
            for (int i = 0; i < 5; i++)
                temp[i] = '-';

        } else if (d <= 5) {
            // 1..5
            int dots   = (int)d;
            int dashes = 5 - dots;

            int i = 0;
            for (; i < dots; i++)
                temp[i] = '.';
            for (; i < 5; i++)
                temp[i] = '-';

        } else { 
            // 6..9
            int dashes = (int)d - 5;

            int i = 0;
            for (; i < dashes; i++)
                temp[i] = '-';
            for (; i < 5; i++)
                temp[i] = '.';
        }

        buff += 5;
        len--;
    }

    return sum;
}
