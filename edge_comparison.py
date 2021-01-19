set_cd = {(2, 1), (4, 3), (6, 5),(34,23),(27,12)}
set_rs = {(1, 2), (3, 4), (5, 6), (7, 8)}
set_true = {(1, 2), (3, 4), (5, 6),(12,27)}

def change_tupel(tuple_list):
    changed_tuple_list = set()
    for x in tuple_list:
        changed_tuple_list.add(((x[1]), (x[0])))
    return changed_tuple_list

# False Positive

tuple_list1 = set_true - set_cd
tuple_list2 = set_true - set_rs

change_tupel1 = change_tupel(tuple_list1)
change_tupel2 = change_tupel(tuple_list2)

result_cd_fp = (set_cd - change_tupel1) - set_true
result_rs_fp = (set_rs - change_tupel2) - set_true

print("result_cd")
print(result_cd_fp)
print("result_rs")
print(result_rs_fp)

## jetzt mal false negative
# vertausche set_true und set_cd/set_rs

tuple_list3 = set_cd - set_true
tuple_list4 = set_rs - set_true

change_tupel3 = change_tupel(tuple_list3)
change_tupel4 = change_tupel(tuple_list4)

result_cd_fn = (set_true - change_tupel3) - set_cd
result_rs_fn = (set_true - change_tupel4) - set_rs

print("result_cd_fn")
print(result_cd_fn)
print("result_rs_fn")
print(result_rs_fn)

#als Letztes true positives
# tp = Vereinigung aus beiden Sets
# Ich verdopple das set_true mit den gedrehten Werten und vereinige beide

set_true_double = set_true | change_tupel(set_true)
tuple_list5 = set_true_double & set_cd
tuple_list6 = set_true_double & set_rs

print(set_true_double)

print("tp cd", tuple_list5)
print("tp rs", tuple_list6)