"""Common functions"""


def process_args(args, ret_type='str', ret_mode=None):
    """Split arguments (only strings) on comma, return in requested

    ret_type
        (str, int or bool)

    as

    ret_mode
        char vector - c(1,2,3)
        list vector - list(1,2,3)
        list - as list [1,2,3]
        None  - same as input

    """
    if args:
        all_args = [item.strip() for item in args.split(',')]
    else:
        all_args = []
    num_options = len(all_args)
    if num_options == 1 and len(all_args[0]):
        if ret_type == 'int':
            option = int(all_args[0])
            if ret_mode == 'charvector':
                # if option:
                return 'c({})'.format(option)
                # else:
                #     return 'c()'
            elif ret_mode == 'listvector':
                # if option:
                return 'list({})'.format(option)
                # else:
                #     return 'list()'
            elif ret_mode == 'list':
                # if option:
                return [option]
                # else:
                #     return []
            else:
                return option
        else:
            # str, bool
            option = all_args[0]
            if ret_mode == 'charvector':
                # if len(option):
                return 'c("{}")'.format(option)
                # else:
                #     return 'c("")'
            elif ret_mode == 'listvector':
                # if option:
                return 'list("{}")'.format(option)
                # else:
                #     return 'list("")'
            elif ret_mode == 'list':
                # if option:
                return [option]
                # else:
                #     return []
            else:
                return option
    elif num_options > 1:
        if ret_type == 'int':
            options = tuple([int(item) for item in all_args])
            if ret_mode == 'charvector':
                # if len(options):
                return 'c{}'.format(options)
                # else:
                #     return 'c()'
            elif ret_mode == 'listvector':
                # if len(options):
                return 'list{}'.format(options)
                # else:
                #     return 'list()'
            elif ret_mode == 'list':
                return list(options)
        else:
            options = tuple(all_args)
            if ret_mode == 'charvector':
                # if len(options):
                return 'c{}'.format(options)
            elif ret_mode == 'listvector':
                return 'list{}'.format(options)
            elif ret_mode == 'list':
                # if len(all_args):
                return all_args
                # else:
                #     return []
            else:
                # as original with spaces stripped
                return ','.join(all_args)
