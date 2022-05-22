#include "ringbuffer.h"


ring_buffer *ring_buffer_new() {
    ring_buffer *my_ring = (ring_buffer *) malloc(sizeof(ring_buffer));
    my_ring->capacity = RING_BUFFER_INIT_SIZE;
    /*use malloc to create a new ring_buffer with length = sizeof(ring_buffer),
     * then assign its capacity with RING_BUFFER_INIT_SIZE.*/
    my_ring->content = (int *) malloc(RING_BUFFER_INIT_SIZE * sizeof(int));
    my_ring->read_pos = my_ring->write_pos = 0;
    return my_ring;
    /*creat the pointer to the content with size equals to RING_BUFFER_INIT_SIZE * sizeof(int) by malloc
     * and set the read_pos and write_pos to the beginning.*/
}

void ring_buffer_destroy(ring_buffer **buffer) {
    if (buffer == NULL || *buffer == NULL) return;
    /*check if the buffer is a null pointer or the *buffer is a null pointer, if they are,
     * then we are done and just return because every memory has been freed.*/
    if ((*buffer)->content != NULL) free((*buffer)->content);
    /*if the content is NULL, we will free content*/
    free(*buffer);
    *buffer = NULL;
    /*if the buffer is not NULL, then firstly we free the content, then we free the buffer
     * and set buffer = NULL to ensure there is no memory leak*/
}

bool ring_buffer_is_empty(const ring_buffer *buffer) {
    if (buffer == NULL) return false;
    if (buffer->read_pos == buffer->write_pos) return true;
    else return false;
    /*According to the piazza, I just regard NULL as not empty and return false
     * then when the buffer is not NULL, we just check whether the read_pos and write_pos are in the same position
     * because when storing those two pointers, I never do "%" operation and the number of those two pointers
     * are exactly the number of the ring_buffer read or wrote historically*/
}

static bool ring_buffer_is_full(const ring_buffer *buffer) {
    if (buffer == NULL) return false;
    if (buffer->write_pos - buffer->read_pos == buffer->capacity) return true;
    else return false;
    /*This is not a function declared in the ringbuffer.h file, but defined by myself, whose purpose is
     * to check whether the ring_buffer is full. Actually the procedure is almost like the ring_buffer_is_empty function
     * and the only difference is when the difference between write_pos and read_pos equals to the capacity
     * can I regard the buffer is full. This is because the difference means how many data the write_pos wrote than the read_pos read
     * if it equals to the capacity, it means the write_pos catches the read_pos up with a circle more*/
}

bool ring_buffer_read(ring_buffer *buffer, int *data) {
    if (buffer == NULL || data == NULL || buffer->content == NULL || ring_buffer_is_empty(buffer)) return false;
    /*If the buffer is NULL or the data is NULL or the buffer is empty content is NULL, it can not read in a new data, so we just return false.*/
    *data = buffer->content[(buffer->read_pos) % buffer->capacity];
    buffer->read_pos++;
    return true;
    /*if the buffer really can read, then we enqueue the data to the position of read_pos % capacity
     * and then add the read_pos with 1 to update the position of read_pos*/
}

bool ring_buffer_write(ring_buffer *buffer, const int data) {
    if (buffer == NULL || buffer->content == NULL) return false;
    /*If the buffer is NULL or the content is NULL, then it can not write, so I just return false.*/
    if (ring_buffer_is_full(buffer)) {
        int *tmp1 = buffer->content, *tmp2;
        size_t i;
        /*if the buffer is full, then I create tmp1 to remember the position the content pointing to
         * and creat tmp2 to be the pointer of new length. Then I create i as iterator*/
        if (buffer->capacity < 1024) {
            tmp2 = (int *) malloc(RING_BUFFER_GROW_FACTOR1 * buffer->capacity * sizeof(int));
        } else {
            tmp2 = (int *) malloc((size_t) (RING_BUFFER_GROW_FACTOR2 * buffer->capacity) * sizeof(int));
        }
        /*if capacity is smaller than 1024, we dilate it with RING_BUFFER_GROW_FACTOR1
         * if the capacity is >= 1024, then we dilate it with RING_BUFFER_GROW_FACTOR2.*/
        for (i = 0; i < buffer->capacity; i++) {
            ring_buffer_read(buffer, &tmp2[i]);
        }
        /*use ring_buffer_read to transit entries from content to tmp2*/
        buffer->read_pos = 0;
        buffer->write_pos = buffer->capacity;
        buffer->content = tmp2;
        /*set the read_pos to 0 because data start in position 0 in the new content
        *update the write_pos and let content equal to tmp2.*/
        free(tmp1);
        if (buffer->capacity < 1024) buffer->capacity *= RING_BUFFER_GROW_FACTOR1;
        else buffer->capacity *= RING_BUFFER_GROW_FACTOR2;
        /*free tmp1 which means free the original content, then update the capacity.*/
    }
    buffer->content[(buffer->write_pos) % buffer->capacity] = data;
    buffer->write_pos++;
    return true;
    /*if the buffer is not full, then enqueue data and update write_pos and return true.*/
}

bool ring_buffer_read_multi(ring_buffer *buffer, size_t rdsize, int *data) {
    size_t i;
    if (buffer == NULL || data == NULL || ((buffer->write_pos - buffer->read_pos) < rdsize) || buffer->content == NULL)
        return false;
    /*Firstly we declare i as the iterator
     * then if buffer is NULL or data is NULL or there is no enough places to store or the content is NULL
     * we just return false without changing anything.*/
    for (i = 0; i < rdsize; ++i) {
        if (!ring_buffer_read(buffer, &data[i])) return false;
    }
    return true;
    /*If we can read in those data, we iterate data and call ring_buffer_read to read in every signal data
     * if the operation return false, then it means there is a problem in reading, and we need to stop the
     * procedure immediately by return false. Else we return true.*/
}

bool ring_buffer_write_multi(ring_buffer *buffer, size_t wrtsize, const int *data) {
    size_t i;
    if (buffer == NULL || data == NULL || buffer->content == NULL) return false;
    /*Firstly we declare i as iterator
     * then if buffer is NULL or data is NULL or the content is NULL, we can not write new data, so we just return false.*/
    for (i = 0; i < wrtsize; ++i) {
        if (!(ring_buffer_write(buffer, data[i]))) return false;
    }
    return true;
    /*Just like the procedure in the ring_buffer_read_multi, we repeatedly call ring_buffer_write to write in every signal data
     * If the operation return false, we need to stop it by return false, else we return true.*/
}

bool ring_buffer_map(ring_buffer *buffer, map_func func) {
    size_t i;
    if (buffer == NULL || func == NULL || ring_buffer_is_empty(buffer) || buffer->content == NULL) return false;
    /*Firstly, we declare i as iterator
     * Then if buffer is NULL or func is NULL or the buffer is empty or the content is NULL, we can not apply any function, so we just return false.*/
    for (i = buffer->read_pos; i < buffer->write_pos; ++i) {
        buffer->content[i % buffer->capacity] = func(buffer->content[i % buffer->capacity]);
    }
    return true;
    /*for every entry in content, we apply the func to it and update its value with the return value*/
}

