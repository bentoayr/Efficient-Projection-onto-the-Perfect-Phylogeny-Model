#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_N_NODES 13
#define topk 100

typedef float realnumber;
typedef unsigned long int longint;
typedef short int shortint;


typedef struct {
    shortint first;
    shortint second;
} edge;

realnumber update_df(realnumber df, realnumber ddf, realnumber t, realnumber t_old) {
    return df + ddf * (t - t_old);
}

void up_to_fixed_ind(shortint start_ind, shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    shortint current_ind = father[start_ind];
    if (current_ind < 0) {
        return;
    }
    *holder_front = 0;
    *holder_back = 0;
    *holder_size = 1;
    holder[*holder_back] = current_ind;
    while (status[current_ind] != 2) {
        current_ind = father[current_ind];
        if (current_ind < 0) {
            return;
        }
        *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
        holder[*holder_back] = current_ind;
        ++(*holder_size);
    }
}

void bfs_to_fixed_get_all_new(shortint start_ind, shortint status[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint node_stack[], shortint *node_stack_back, shortint *node_stack_size) {
    
    *queue_front = 0;
    queue[*queue_front] = start_ind;
    *queue_back = 0;
    *queue_size = 1;
    while (*queue_size != 0) {
        shortint current_ind = queue[*queue_front];
        *queue_front = (*queue_front + 1) % (MAX_N_NODES + 1);
        --*queue_size;
        *node_stack_back = (*node_stack_back + 1) % (MAX_N_NODES + 1);
        node_stack[*node_stack_back] = current_ind;
        ++(*node_stack_size);
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            if (status[kid_ind] != 2 && status[kid_ind] != 3) {
                *queue_back = (*queue_back + 1) % (MAX_N_NODES + 1);
                queue[*queue_back] = kid_ind;
                ++*queue_size;
            }
        }
    }
}


realnumber update_ddf_reborn(realnumber ddf, shortint new_fixed_ind, realnumber v[], realnumber v_old[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim, shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size, shortint seen[], realnumber gamma[]) {
    realnumber old_contribution = 0;
    realnumber new_contribution = 0;
    
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    shortint current_root_ind = new_fixed_ind;
    if (father[new_fixed_ind] == 0) {
        old_contribution += (v_old[new_fixed_ind] - 0) * (v_old[new_fixed_ind] - 0);
        new_contribution += (v[new_fixed_ind] - 0) * (v[new_fixed_ind] - 0);
    }
    else {
        up_to_fixed_ind(new_fixed_ind, status, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
        current_root_ind = holder[(*holder_back - 1) % (MAX_N_NODES + 1)];
        
        old_contribution += (v_old[current_root_ind] - v_old[father[current_root_ind]]) * (v_old[current_root_ind] - v_old[father[current_root_ind]]);
        new_contribution += (v[current_root_ind] - v[father[current_root_ind]]) * (v[current_root_ind] - v[father[current_root_ind]]);
        
        *holder_front = 0;
        *holder_back = (MAX_N_NODES + 1) - 1;
        *holder_size = 0;
        
        bfs_to_fixed_get_all_new(current_root_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    }
    
    shortint current_num_kids = num_kids[new_fixed_ind];
    for (shortint i = 0; i < current_num_kids; ++i) {
        shortint current_kid_ind = kids[new_fixed_ind * (MAX_N_NODES + 1) + i];
        bfs_to_fixed_get_all_new(current_kid_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
        
        old_contribution += (v_old[current_kid_ind] - v_old[father[current_kid_ind]]) * (v_old[current_kid_ind] - v_old[father[current_kid_ind]]);
        new_contribution += (v[current_kid_ind] - v[father[current_kid_ind]]) * (v[current_kid_ind] - v[father[current_kid_ind]]);
    }
    
    while (*holder_size != 0) {
        shortint current_node_ind = holder[*holder_front];
        *holder_front = (*holder_front + 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        
        shortint kids_num = num_kids[current_node_ind];
        for (shortint i = 0; i < kids_num; ++i) {
            old_contribution += (v_old[current_node_ind] - v_old[kids[current_node_ind * (MAX_N_NODES + 1) + i]]) * (v_old[current_node_ind] - v_old[kids[current_node_ind * (MAX_N_NODES + 1) + i]]);
            new_contribution += (v[current_node_ind] - v[kids[current_node_ind * (MAX_N_NODES + 1) + i]]) * (v[current_node_ind] - v[kids[current_node_ind * (MAX_N_NODES + 1) + i]]);
        }
        
        seen[current_node_ind] = 0;
        
        status[current_node_ind] = 0;
        
        gamma[current_node_ind] = 1;
    }
    
    return ddf + (new_contribution - old_contribution);
}

void heap_insert(realnumber intersetion, shortint node_index, realnumber max_heap[], shortint heap_to_tree_index[], shortint tree_to_heap_index[], shortint *heap_size) {
    ++(*heap_size);
    max_heap[*heap_size - 1] = intersetion;
    heap_to_tree_index[*heap_size - 1] = node_index;
    tree_to_heap_index[node_index] = *heap_size - 1;
    shortint current_ind = *heap_size - 1;
    while (1) {
        shortint father_ind = (current_ind - 1) / 2;
        if (intersetion > max_heap[father_ind]) {
            max_heap[current_ind] = max_heap[father_ind];
            max_heap[father_ind] = intersetion;
            heap_to_tree_index[current_ind] = heap_to_tree_index[father_ind];
            heap_to_tree_index[father_ind] = node_index;
            tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
            tree_to_heap_index[node_index] = father_ind;
            current_ind = father_ind;
        }
        else {
            break;
        }
    }
}





void heap_delete(shortint node_index, realnumber max_heap[], shortint heap_to_tree_index[], shortint tree_to_heap_index[], shortint *heap_size) {
    
    
    shortint current_ind = tree_to_heap_index[node_index];
    realnumber temp = max_heap[*heap_size - 1];
    max_heap[current_ind] = temp;
    shortint temp_node_ind = heap_to_tree_index[*heap_size - 1];
    heap_to_tree_index[current_ind] = temp_node_ind;
    tree_to_heap_index[temp_node_ind] = current_ind;
    
    --(*heap_size);
    if (max_heap[current_ind] > max_heap[(current_ind - 1) / 2]) {
        
        while (1) {
            shortint father_ind = (current_ind - 1) / 2; // Integer division acts like floor function
            if (max_heap[current_ind] > max_heap[father_ind]) {
                max_heap[current_ind] = max_heap[father_ind];
                max_heap[father_ind] = temp;
                heap_to_tree_index[current_ind] = heap_to_tree_index[father_ind];
                heap_to_tree_index[father_ind] = temp_node_ind;
                tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                tree_to_heap_index[temp_node_ind] = father_ind;
                current_ind = father_ind;
            }
            else {
                break;
            }
        }
    }
    else {
        
        while (1) {
            shortint left_ind = 2 * current_ind + 1;
            shortint right_ind = 2 * current_ind + 2;
            
            if (left_ind < *heap_size && right_ind < *heap_size) {
                if (temp < max_heap[left_ind] && temp < max_heap[right_ind]) {
                    if (max_heap[left_ind] >= max_heap[right_ind]) {
                        max_heap[current_ind] = max_heap[left_ind];
                        max_heap[left_ind] = temp;
                        heap_to_tree_index[current_ind] = heap_to_tree_index[left_ind];
                        heap_to_tree_index[left_ind] = temp_node_ind;
                        tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                        tree_to_heap_index[temp_node_ind] = left_ind;
                        current_ind = left_ind;
                    }
                    else {
                        max_heap[current_ind] = max_heap[right_ind];
                        max_heap[right_ind] = temp;
                        heap_to_tree_index[current_ind] = heap_to_tree_index[right_ind];
                        heap_to_tree_index[right_ind] = temp_node_ind;
                        tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                        tree_to_heap_index[temp_node_ind] = right_ind;
                        current_ind = right_ind;
                    }
                }
                else if (temp < max_heap[left_ind] && temp >= max_heap[right_ind]) {
                    max_heap[current_ind] = max_heap[left_ind];
                    max_heap[left_ind] = temp;
                    heap_to_tree_index[current_ind] = heap_to_tree_index[left_ind];
                    heap_to_tree_index[left_ind] = temp_node_ind;
                    tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                    tree_to_heap_index[temp_node_ind] = left_ind;
                    current_ind = left_ind;
                }
                else if (temp < max_heap[right_ind] && temp >= max_heap[left_ind]) {
                    max_heap[current_ind] = max_heap[right_ind];
                    max_heap[right_ind] = temp;
                    heap_to_tree_index[current_ind] = heap_to_tree_index[right_ind];
                    heap_to_tree_index[right_ind] = temp_node_ind;
                    tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                    tree_to_heap_index[temp_node_ind] = right_ind;
                    current_ind = right_ind;
                }
                else {
                    break;
                }
            }
            else if (left_ind < *heap_size && temp < max_heap[left_ind]) {
                max_heap[current_ind] = max_heap[left_ind];
                max_heap[left_ind] = temp;
                heap_to_tree_index[current_ind] = heap_to_tree_index[left_ind];
                heap_to_tree_index[left_ind] = temp_node_ind;
                tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                tree_to_heap_index[temp_node_ind] = left_ind;
                current_ind = left_ind;
            }
            else if (right_ind < *heap_size && temp < max_heap[right_ind]) {
                max_heap[current_ind] = max_heap[right_ind];
                max_heap[right_ind] = temp;
                heap_to_tree_index[current_ind] = heap_to_tree_index[right_ind];
                heap_to_tree_index[right_ind] = temp_node_ind;
                tree_to_heap_index[heap_to_tree_index[current_ind]] = current_ind;
                tree_to_heap_index[temp_node_ind] = right_ind;
                current_ind = right_ind;
            }
            else {
                break;
            }
        }
    }
    
}


shortint next_turn_better_new(shortint new_fixed_ind, realnumber t_pre, realnumber z[], realnumber v[], realnumber value[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim, shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size, realnumber max_heap[], shortint heap_to_tree_index[], shortint tree_to_heap_index[], shortint *heap_size) {
    v[0] = 0;
    
    if (tree_to_heap_index[new_fixed_ind] < *heap_size) {
        heap_delete(new_fixed_ind, max_heap, heap_to_tree_index, tree_to_heap_index, heap_size);
    }
    
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    if (father[new_fixed_ind] != 0) {
        up_to_fixed_ind(new_fixed_ind, status, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
        shortint current_root_ind = holder[(*holder_back - 1) % (MAX_N_NODES + 1)];
        
        heap_delete(current_root_ind, max_heap, heap_to_tree_index, tree_to_heap_index, heap_size);
        
        *holder_front = 0;
        *holder_back = (MAX_N_NODES + 1) - 1;
        *holder_size = 0;
        
        bfs_to_fixed_get_all_new(current_root_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
        realnumber t_possilbe = -INFINITY;
        while (*holder_size != 0) {
            shortint current_node_ind = holder[*holder_front];
            *holder_front = (*holder_front + 1) % (MAX_N_NODES + 1);
            --(*holder_size);
            
            if (status[current_node_ind] != 2) {
                realnumber curr_t_possilbe = (z[current_node_ind] - v[current_node_ind] * t_pre + value[current_node_ind]) / (1 - v[current_node_ind]);
                if (t_possilbe < curr_t_possilbe) {
                    t_possilbe = curr_t_possilbe;
                }
            }
        }
        
        heap_insert(t_possilbe, current_root_ind, max_heap, heap_to_tree_index, tree_to_heap_index, heap_size);
    }
    
    int current_num_kids = num_kids[new_fixed_ind];
    for (shortint i = 0; i < current_num_kids; ++i) {
        
        *holder_front = 0;
        *holder_back = (MAX_N_NODES + 1) - 1;
        *holder_size = 0;
        shortint current_kid_ind = kids[new_fixed_ind * (MAX_N_NODES + 1) + i];
        realnumber t_possilbe = -INFINITY;
        if (status[current_kid_ind] != 2) {
            bfs_to_fixed_get_all_new(current_kid_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
            while (*holder_size != 0) {
                shortint current_node_ind = holder[*holder_front];
                *holder_front = (*holder_front + 1) % (MAX_N_NODES + 1);
                --(*holder_size);
                
                if (status[current_node_ind] != 2) {
                    realnumber curr_t_possilbe = (z[current_node_ind] - v[current_node_ind] * t_pre + value[current_node_ind]) / (1 - v[current_node_ind]);
                    if (t_possilbe < curr_t_possilbe) {
                        t_possilbe = curr_t_possilbe;
                    }
                }
            }
        }
        
        heap_insert(t_possilbe, current_kid_ind, max_heap, heap_to_tree_index, tree_to_heap_index, heap_size);
    }
    
    shortint max_subtree_root = heap_to_tree_index[0];
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    bfs_to_fixed_get_all_new(max_subtree_root, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    realnumber t_want = max_heap[0];
    while (*holder_size != 0) {
        shortint current_node_ind = holder[*holder_front];
        *holder_front = (*holder_front + 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        // Find max intersection:
        if (status[current_node_ind] != 2) {
            realnumber curr_t_possilbe = (z[current_node_ind] - v[current_node_ind] * t_pre + value[current_node_ind]) / (1 - v[current_node_ind]);
            if (t_want == curr_t_possilbe) {
                new_fixed_ind = current_node_ind;
            }
        }
    }
    return new_fixed_ind;
}

void prune(shortint root, shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    
    *queue_front = 0;
    queue[*queue_front] = root;
    *queue_back = 0;
    *queue_size = 1;
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    shortint leaf_count = 0;
    
    while (*queue_size != 0) {
        shortint current_ind = queue[*queue_front];
        *queue_front = (*queue_front + 1) % (MAX_N_NODES + 1);
        --*queue_size;
        *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
        holder[*holder_back] = current_ind;
        ++(*holder_size);
        if (num_kids[current_ind] == 0) {
            ++leaf_count;
        }
        else if (num_kids[current_ind] == 1 && status[current_ind] == 2) {
            ++leaf_count;
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            if (status[kid_ind] != 2) {
                *queue_back = (*queue_back + 1) % (MAX_N_NODES + 1);
                queue[*queue_back] = kid_ind;
                ++*queue_size;
            }
        }
    }
    
    while (*holder_size != 0) {
        
        shortint current_ind = holder[*holder_back];
        *holder_back = (*holder_back - 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        
        if (status[current_ind] != 2) {
            
            if (leaf_count != 0) {
                --leaf_count;
                if (num_kids[current_ind] == 0) {
                    status[current_ind] = 1;
                }
            }
            
            else {
                shortint unstretched_indicator = 1;
                for (shortint i = 0; i < num_kids[current_ind]; ++i) {
                    shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
                    if (status[kid_ind] != 1) {
                        unstretched_indicator = 0;
                        break;
                    }
                }
                if (unstretched_indicator == 1) {
                    status[current_ind] = 1;
                }
            }
        }
    }
}



void update_v_reborn_stretched(shortint root, realnumber v[], realnumber gamma[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint seen[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    
    *queue_front = 0;
    queue[*queue_front] = root;
    *queue_back = 0;
    *queue_size = 1;
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    while (*queue_size != 0) {
        shortint current_ind = queue[*queue_front];
        *queue_front = (*queue_front + 1) % (MAX_N_NODES + 1);
        --(*queue_size);
        
        if (status[current_ind] == 0) {
            *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
            holder[*holder_back] = current_ind;
            ++(*holder_size);
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            if (status[kid_ind] != 2) {
                *queue_back = (*queue_back + 1) % (MAX_N_NODES + 1);
                queue[*queue_back] = kid_ind;
                ++*queue_size;
            }
        }
    }
    
    while (*holder_size != 1) {
        int current_ind = holder[*holder_back];
        *holder_back = (*holder_back - 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        seen[current_ind] = 1;
        realnumber numerator = 0, denominator = 0;
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            
            if (status[kid_ind] != 1) {
                numerator += gamma[kid_ind] * v[kid_ind];
                denominator += gamma[kid_ind];
            }
        }
        v[current_ind] = numerator / denominator;
        gamma[current_ind] = 1 / ((1 / gamma[current_ind]) + (1 / denominator));
    }
    
    *holder_front = 0;
    *holder_back = 0;
    holder[*holder_back] = root;
    *holder_size = 1;
    
    while (*holder_size != 0) {
        shortint current_ind = holder[*holder_back];
        *holder_back = (*holder_back - 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        
        shortint parent_ind = father[current_ind];
        gamma[current_ind] = 1;
        realnumber numerator = gamma[current_ind] * v[parent_ind], denominator = gamma[current_ind];
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            
            if (status[kid_ind] != 1) {
                numerator += gamma[kid_ind] * v[kid_ind];
                denominator += gamma[kid_ind];
            }
        }
        v[current_ind] = numerator / denominator;
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            
            if (status[kid_ind] == 0) {
                *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
                holder[*holder_back] = kid_ind;
                ++(*holder_size);
            }
        }
    }
}

void update_v_reborn_unstretched(shortint root, realnumber v[], realnumber gamma[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size) {
    
    *queue_front = 0;
    queue[*queue_front] = root;
    *queue_back = 0;
    *queue_size = 1;
    
    while (*queue_size != 0) {
        shortint current_ind = queue[*queue_front];
        *queue_front = (*queue_front + 1) % (MAX_N_NODES + 1);
        --*queue_size;
        if (status[current_ind] == 1) {
            v[current_ind] = v[father[current_ind]];
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            if (status[kid_ind] != 2) {
                *queue_back = (*queue_back + 1) % (MAX_N_NODES + 1);
                queue[*queue_back] = kid_ind;
                ++*queue_size;
            }
        }
    }
}


void update_v_reborn(shortint new_fixed_ind, realnumber v[], realnumber gamma[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size, shortint seen[]) {
    v[0] = 0;
    v[new_fixed_ind] = 1;
    
    *queue_front = 0;
    *queue_back = (MAX_N_NODES + 1) - 1;
    *queue_size = 0;
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    if (father[new_fixed_ind] != 0) {
        up_to_fixed_ind(new_fixed_ind, status, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
        shortint current_root_ind = holder[(*holder_back - 1) % (MAX_N_NODES + 1)];
        
        prune(current_root_ind, status, father, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_front, holder_back, holder_size);
        
        update_v_reborn_stretched(current_root_ind, v, gamma, status, father, num_kids, kids, seen, queue, queue_front, queue_back, queue_size, holder, holder_front, holder_back, holder_size);
        
        update_v_reborn_unstretched(current_root_ind, v, gamma, status, father, num_kids, kids, queue, queue_front, queue_back, queue_size);
    }
    
    for (shortint i = 0; i < num_kids[new_fixed_ind]; ++i) {
        shortint current_root_ind = kids[new_fixed_ind * (MAX_N_NODES + 1) + i];
        
        prune(current_root_ind, status, father, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_front, holder_back, holder_size);
        
        if (status[current_root_ind] == 0) {
            update_v_reborn_stretched(current_root_ind, v, gamma, status, father, num_kids, kids, seen, queue, queue_front, queue_back, queue_size, holder, holder_front, holder_back, holder_size);
        }
        
        update_v_reborn_unstretched(current_root_ind, v, gamma, status, father, num_kids, kids, queue, queue_front, queue_back, queue_size);
    }
    v[0] = 0;
    v[new_fixed_ind] = 1;
}


void bfs_get_level_new(shortint status[], shortint level[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size) {
    
    *queue_front = 0;
    queue[*queue_front] = 0;
    *queue_back = 0;
    *queue_size = 1;
    level[0] = 0;
    
    while (*queue_size != 0) {
        shortint current_ind = queue[*queue_front];
        *queue_front = (*queue_front + 1) % (MAX_N_NODES + 1);
        --*queue_size;
        
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            *queue_back = (*queue_back + 1) % (MAX_N_NODES + 1);
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            queue[*queue_back] = kid_ind;
            ++*queue_size;
            level[kid_ind] = level[current_ind] + 1;
        }
    }
}

realnumber dfs_tree_cost_from_z_non_recursive(shortint num_nodes, shortint *adj_list, shortint *final_degrees, shortint root_node, shortint *stack, shortint *visited, realnumber *z, realnumber *data, shortint T, shortint t) {
    
    shortint curr = root_node;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    
    realnumber temp2 = (z[curr + 1] + data[curr*T + t]);
    realnumber temp = (temp2*temp2);
    
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < final_degrees[curr]; j++){
            shortint child_ix = adj_list[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root_node]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                temp2 = ((z[child_ix + 1] - z[curr + 1]) + data[child_ix*T + t]);
                temp +=  (temp2*temp2);
            }
        }
    }
    return temp;
}

void update_z_new(shortint new_fixed_ind, realnumber value[], realnumber t, realnumber z[], realnumber z_timestamp[], realnumber v[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    if (father[new_fixed_ind] != 0) {
        up_to_fixed_ind(new_fixed_ind, status, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
        shortint current_root_ind = holder[(*holder_back - 1) % (MAX_N_NODES + 1)]; // The second last element of holder is the root of the current subtree
        
        *holder_front = 0;
        *holder_back = (MAX_N_NODES + 1) - 1;
        *holder_size = 0;
        bfs_to_fixed_get_all_new(current_root_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    }
    
    shortint current_num_kids = num_kids[new_fixed_ind];
    for (shortint i = 0; i < current_num_kids; ++i) {
        shortint current_root_ind = kids[new_fixed_ind * (MAX_N_NODES + 1) + i];
        bfs_to_fixed_get_all_new(current_root_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    }
    
    while (*holder_size != 0) {
        shortint ind = holder[*holder_back];
        *holder_back = (*holder_back - 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        z[ind] = z[ind] + v[ind] * (t - z_timestamp[ind]);
        z_timestamp[ind] = t;
    }
    z[new_fixed_ind] = z[new_fixed_ind] + v[new_fixed_ind] * (t - z_timestamp[new_fixed_ind]);
    z_timestamp[new_fixed_ind] = t;
}

void update_v_old_new(shortint new_fixed_ind, shortint pre_fixed_ind, realnumber t, realnumber v[], realnumber v_old[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim, shortint queue[], shortint *queue_front, shortint *queue_back, shortint *queue_size, shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    
    *holder_front = 0;
    *holder_back = (MAX_N_NODES + 1) - 1;
    *holder_size = 0;
    
    if (father[new_fixed_ind] != 0) {
        up_to_fixed_ind(new_fixed_ind, status, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
        shortint current_root_ind = holder[(*holder_back - 1) % (MAX_N_NODES + 1)];
        
        *holder_front = 0;
        *holder_back = (MAX_N_NODES + 1) - 1;
        *holder_size = 0;
        
        bfs_to_fixed_get_all_new(current_root_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    }
    
    shortint current_num_kids = num_kids[new_fixed_ind];
    for (shortint i = 0; i < current_num_kids; ++i) {
        shortint current_kid_ind = kids[new_fixed_ind * (MAX_N_NODES + 1) + i];
        bfs_to_fixed_get_all_new(current_kid_ind, status, num_kids, kids, queue, queue_front, queue_back, queue_size, holder, holder_back, holder_size);
    }
    
    while (*holder_size != 0) {
        shortint current_node_ind = holder[*holder_front];
        *holder_front = (*holder_front + 1) % (MAX_N_NODES + 1);
        --(*holder_size);
        v_old[current_node_ind] = v[current_node_ind];
    }
    v_old[new_fixed_ind] = v[new_fixed_ind];
    v_old[pre_fixed_ind] = v[pre_fixed_ind];
}


void best_tree_reborn(realnumber z[], realnumber value[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim) {
    
    shortint queue[(MAX_N_NODES + 1)];
    shortint queue_front = 0;
    shortint queue_back = (MAX_N_NODES + 1) - 1;
    shortint queue_size = 0;
    shortint holder[(MAX_N_NODES + 1)];
    shortint holder_front = 0;
    shortint holder_back = (MAX_N_NODES + 1) - 1;
    shortint holder_size = 0;
    
    shortint level[(MAX_N_NODES + 1)];
    bfs_get_level_new(status, level, num_kids, kids, queue, &queue_front, &queue_back, &queue_size);
    
    shortint seen[(MAX_N_NODES + 1)] = {0};
    realnumber v[(MAX_N_NODES + 1)] = {0};
    realnumber v_old[(MAX_N_NODES + 1)] = {0};
    realnumber gamma[(MAX_N_NODES + 1)] = {1};
    realnumber max_heap[(MAX_N_NODES + 1)];
    shortint heap_to_tree_index[(MAX_N_NODES + 1)];
    shortint tree_to_heap_index[(MAX_N_NODES + 1)];
    for (shortint i = 0; i <= dim; ++i) {
        max_heap[i] = -INFINITY;
        heap_to_tree_index[i] = kids[0];
        tree_to_heap_index[i] = kids[0];
        gamma[i] = 1;
        z[i] = 0;
        v[i] = 0;
    }
    shortint heap_size = 0;
    
    realnumber n_max = value[1];
    shortint new_fixed_ind = 0, pre_fixed_ind = 0;
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] >= n_max) {
            n_max = value[i];
        }
    }
    shortint top_level = dim + 1;
    
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] == n_max && level[i] < top_level) {
            new_fixed_ind = i;
            top_level = level[i];
        }
    }
    
    realnumber t = n_max, t_old = n_max;
    heap_insert(t, kids[0], max_heap, heap_to_tree_index, tree_to_heap_index, &heap_size);
    
    status[new_fixed_ind] = 2;
    seen[new_fixed_ind] = 1;
    update_v_reborn(new_fixed_ind, v, gamma, status, father, num_kids, kids, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size, seen);
    realnumber z_timestamp[MAX_N_NODES + 1];
    
    for (shortint i = 0; i <= dim; ++i) {
        z_timestamp[i] = t;
        seen[i] = 0;
    }
    
    seen[new_fixed_ind] = 1;
    z[new_fixed_ind] = z[new_fixed_ind] + (t - z_timestamp[new_fixed_ind]);
    z_timestamp[new_fixed_ind] = t;
    realnumber df = 0, df_old = 0, ddf = 0;
    
    while (df > -1) {
        t_old = t;
        df_old = df;
        pre_fixed_ind = new_fixed_ind;
        
        new_fixed_ind = next_turn_better_new(new_fixed_ind, t_old, z, v, value, status, father, num_kids, kids, dim, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size, max_heap, heap_to_tree_index, tree_to_heap_index, &heap_size);
        t = max_heap[0];
        if (t == -INFINITY) {
            ddf = update_ddf_reborn(ddf, pre_fixed_ind, v, v_old, status, father, num_kids, kids, dim, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size, seen, gamma);
            break;
        }
        
        ddf = update_ddf_reborn(ddf, pre_fixed_ind, v, v_old, status, father, num_kids, kids, dim, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size, seen, gamma);
        
        df = update_df(df, ddf, t, t_old);
        
        status[new_fixed_ind] = 2;
        seen[new_fixed_ind] = 1;
        
        if (df > -1) {
            update_z_new(new_fixed_ind, value, t, z, z_timestamp, v, status, father, num_kids, kids, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size);
            
            update_v_old_new(new_fixed_ind, pre_fixed_ind, t, v, v_old, status, father, num_kids, kids, dim, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size);
            
            update_v_reborn(new_fixed_ind, v, gamma, status, father, num_kids, kids, queue, &queue_front, &queue_back, &queue_size, holder, &holder_front, &holder_back, &holder_size, seen);
        }
    }
    
    t = t_old - (df_old + 1) / ddf;
    z[0] = t;
    for (shortint i = 1; i <= dim; ++i) {
        z[i] = z[i] + v[i] * (t - z_timestamp[i]);
    }
}


void dfs_tree_compute_ntilde_non_recursive_array(shortint num_nodes, shortint *adj, shortint *deg, shortint root, shortint *stack, shortint *visited, realnumber *ntilde, realnumber *data,  shortint T, shortint t){
    ntilde[0] = 0;
    shortint curr = root;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    ntilde[curr + 1] = data[curr*T + t];
    
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < deg[curr]; j++){
            shortint child_ix = adj[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                ntilde[child_ix + 1] = ntilde[curr + 1] + data[child_ix*T + t];
            }
        }
    }
}


void convert_tree_data(shortint num_nodes, shortint *final_degrees, shortint *adj_list, shortint *father_list, shortint root_node, shortint num_kids[], shortint kids[]) {
    num_kids[0] = 1;
    kids[0] = root_node + 1;
    for (int i = 1; i <= num_nodes; ++i) {
        if (i == root_node + 1) {
            num_kids[i] = final_degrees[i - 1];
        }
        else {
            num_kids[i] = final_degrees[i - 1] - 1;
        }
    }
    for (shortint i = 1; i <= num_nodes; ++i) {
        if (i == root_node + 1) {
            for (shortint j = 0; j <= num_kids[i]; ++j) {
                if (j == num_kids[i]) {
                    kids[i * (MAX_N_NODES + 1) + num_kids[i]] = -2;
                }
                else {
                    kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j] + 1;
                }
            }
        }
        else {
            shortint father_ind = father_list[i];
            shortint j = 0;
            while (j < num_kids[i] && adj_list[(i-1)*num_nodes + j] + 1 != father_ind) {
                kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j] + 1;
                ++j;
            }
            while (j < num_kids[i]) {
                kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j + 1] + 1;
                ++j;
            }
            if (j == num_kids[i]) {
                kids[i * (MAX_N_NODES + 1) + num_kids[i]] = -2;
            }
        }
    }
    
}

void dfs_tree_compute_fathers_non_recursive_array(shortint num_nodes,shortint *fathers_list, shortint *adj, shortint *deg, shortint root, shortint *stack, shortint *visited){
    fathers_list[0] = -1;
    shortint curr = root;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    fathers_list[root + 1] = 0;
    
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < deg[curr]; j++){
            shortint child_ix = adj[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                fathers_list[child_ix + 1] = curr + 1;
            }
        }
    }
}


realnumber tree_cost_projection(shortint num_nodes, shortint T, realnumber *data, shortint root_node, edge *tree, shortint *adjacency_mat, shortint *final_degrees, shortint *adj_list){
    
    realnumber z[2*MAX_N_NODES + 1];
    realnumber ntilde[2*MAX_N_NODES + 1];
    shortint stack[MAX_N_NODES];
    shortint visited[MAX_N_NODES];
    shortint father_list[2*MAX_N_NODES + 1];
    shortint num_kids[MAX_N_NODES + 1];
    shortint kids[(MAX_N_NODES + 1) * (MAX_N_NODES + 1)];
    
    for (shortint i = 0; i < num_nodes; i++){
        visited[i] = 0;
    }
    
    realnumber error_tree_model = 0.0;
    
    dfs_tree_compute_fathers_non_recursive_array(num_nodes, father_list, adj_list, final_degrees, root_node, stack, visited);
    convert_tree_data(num_nodes, final_degrees, adj_list, father_list, root_node, num_kids, kids);
    
    for (shortint t = 0; t < T; ++t) {
        shortint status[2*MAX_N_NODES + 1] = {0};
        status[0] = 2;
        
        dfs_tree_compute_ntilde_non_recursive_array(num_nodes, adj_list, final_degrees, root_node, stack, visited, ntilde, data, T, t);
        
        best_tree_reborn(z, ntilde, status, father_list, num_kids, kids, num_nodes);
        
        realnumber tmp_cost = dfs_tree_cost_from_z_non_recursive(num_nodes, adj_list, final_degrees, root_node, stack, visited, z, data, T, t);
        
        error_tree_model += tmp_cost;
    }
    
    return error_tree_model;
}


int main(int argc, const char * argv[]) {
    
    // number of mutations, q in the paper
    shortint num_nodes = 10;
    // number of samples, p in the paper.
    shortint T = 1;
    // this is the input, F_hat in the paper with
    realnumber data[10] = {0.11,0.42,0.99,0.84,0.15,0.01,0.14,0.55,0.28,0.46};
    // this is the root of the tree
    shortint root_node = 2;
    
    // list of edges in the tree
    edge tree_edges[9] = {{2,0},{0,1},{1,5},{0,6},{0,7},{7,8},{9,8},{8,3},{3,4}};
    // adjacency matrixx of the tree
    shortint adjacency_mat[10*10] = {0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,
        0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
        0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,
        0,0,0,0,0,0,1,0};
    // list of degrees of each node in the tree
    shortint final_degrees[10] = {4, 2, 1, 2, 1, 1, 1, 2, 3, 1};
    // adjacency list representation of the tree. Each row is the list of neigbours of a node.
    
    shortint adj_list[10*10] = {2, 1, 7, 6,-1,-1,-1,-1,-1,-1
        , 5, 0,-1,-1,-1,-1,-1,-1,-1,-1
        , 0,-1,-1,-1,-1,-1,-1,-1,-1,-1
        , 4, 8,-1,-1,-1,-1,-1,-1,-1,-1
        , 3,-1,-1,-1,-1,-1,-1,-1,-1,-1
        , 1,-1,-1,-1,-1,-1,-1,-1,-1,-1
        , 0,-1,-1,-1,-1,-1,-1,-1,-1,-1
        , 0, 8,-1,-1,-1,-1,-1,-1,-1,-1
        , 7, 3, 9,-1,-1,-1,-1,-1,-1,-1
        , 8,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    
    realnumber cost = tree_cost_projection(num_nodes, T, data,  root_node, tree_edges, adjacency_mat, final_degrees, adj_list);
    
    printf("Cost = %f\n", cost);
    
    return 0;
}

