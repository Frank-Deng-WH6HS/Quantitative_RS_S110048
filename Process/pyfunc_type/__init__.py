#函数类型强制检查服务

import inspect, types, typing; 

#逐一比对传入函数中实参类型与声明的形参类型
def _param_type_enum(func, params, arg_config): 
    for param in params.items(): 
        param_name = param[0]; 
        param_type_decl = param[1].annotation; 
        #如果某个参数没有在annotation说明类型, 则不检查该参数
        if param_type_decl is inspect._empty: 
            continue; 
        #获取函数调用时传入的实参类型
        param_val_passed = arg_config[param[0]]; 
        param_type_passed = type(param_val_passed); 
        #参数类型不相符时报错, 说明正确的类型和当前传入的错误类型
        if type(param_type_decl) is type: 
            if param_type_passed is not param_type_decl: 
                raise TypeError("%s: '%s' must be %s, not %s" % (
                        func.__name__, 
                        param_name, 
                        param_type_decl, 
                        param_type_passed
                    ) ); 
        if type(param_type_decl) is set: 
            if param_type_passed not in param_type_decl: 
                raise TypeError("%s: '%s' must be %s, not %s" % (
                        func.__name__, 
                        param_name, 
                        "\x20or\x20".join(
                            str(elem) for elem in param_type_decl
                        ), 
                        param_type_passed
                    ) ); 

#检查函数中传入参数的类型
#用法: 在函数的def语句之后使用
def type_check(): 
    #读取当前type_check()被调用时所在的函数和模块名称
    context = inspect.getouterframes(inspect.currentframe()); 
    func_context = context[1].frame.f_globals["__name__"]; 
    import __main__; 
    if func_context == __main__.__name__: 
        func = context[1].function; 
        func = context[2].frame.f_globals[func]; 
    else: 
        func = ".".join([__main__.__name__, 
            func_context, context[1].function]); 
        func = eval(func); 
    #读取函数声明时的形式参数名称和annotation类型
    params = inspect.signature(func).parameters; 
    #读取函数调用时的传入的实参值
    arg_config = context[1].frame.f_locals; 
    #逐一检查类型
    _param_type_enum(func, params, arg_config)
    