struct Node {
    long data;
    struct Node *next;
};

/* NOTE: (void *) 0 is how NULL is defined in C.
   Basically, NULL should correspond to zero in your
   y86-64 code, nothing confusing! */

/* Merge the given two sorted linked lists into a single sorted list using
 * recursive Merge operation.
 * Return NULL if both given lists are empty. */
struct Node* merge_sorted_rec(struct Node *a, struct Node *b)
{
    if (a == NULL)
        return b;
    if (b == NULL)
        return a;

    if (a->data <= b->data) {
        a->next = merge_sorted_rec(a->next, b);
        return a;
    } else {
        b->next = merge_sorted_rec(a, b->next);
        return b;
    }
}

/* Merge the two given sorted linked lists into a single sorted list
 * using iterative Merge operation.
 * Return NULL if both linked lists are empty. */
struct Node* merge_sorted_it(struct Node *a, struct Node *b)
{
    struct Node dummy;
    struct Node *tail = &dummy;
    dummy.next = NULL;

    while (a != NULL && b != NULL) {
        if (a->data <= b->data) {
            tail->next = a;
            a = a->next;
        } else {
            tail->next = b;
            b = b->next;
        }
        tail = tail->next;
    }

    if (a != NULL)
        tail->next = a;
    else
        tail->next = b;

    return dummy.next;
}
