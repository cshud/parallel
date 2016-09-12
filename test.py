# import mpi4py.MPI as MPI
#
#
# def pa(a, b, n):
#         # if comm_rank == 0:
#         # temp = range(20)
#         # k = len(temp) / comm_size
#         # data = [temp[i * k:i * k + k] for i in range(comm_size)]
#         # print data
#         #
#         # local_data = comm.scatter(data, root=0)
#         # result = [2*i for i in local_data]
#         # print 'rank %d, got and do:' % comm_rank
#         # print local_data
#         # combine_data = comm.gather(result, root=0)
#         # if comm_rank == 0:
#         # print combine_data
#     x = [0]*n
#     # local_a = comm.scatter([a] * comm_size, root=0)
#     # local_b = comm.scatter([b] * comm_size, root=0)
#
#     for i in range(n):
#         x[i] = a +b*2
#
#
#     combine_a = comm.gather(x, root=0)
#     if comm_rank == 0:
#         print combine_a
#
# #
# comm = MPI.COMM_WORLD
# comm_rank = comm.Get_rank()
# comm_size = comm.Get_size()
# # # if comm_rank == 0:
# # # a = range(100)
# # #     size = 4
# # #     l = len(a)/size
# # #     result = [a[i*l:i*l+l] for i in range(size)]
# # #     print result
# # # pa()
# pa(1, 2, 5)
# # a=1
# # b=2
# # n=5
# # x=[0]*n
# # for it in range(10):
# #     for i in range(n):
# #         x[i] = a+b*2
# # print x