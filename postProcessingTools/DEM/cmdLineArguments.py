import os
import sys

class cmdLineArgs:
    """This class interprets the command-line arguments for post-processing""" + """programs."""

    def __init__(self):
        self.attributes = {
            "root_paths" : [],
            "labels" : {},
            "titles" : [],
            "linestyles" : {},
            "linewidths" : {},
            "colors" : {}
        }
        self.attributeSymbols = {
            "-l" : "labels",
            "-t" : "titles",
            "-ls" : "linestyles",
            "-lw" : "linewidths",
            "-c" : "colors"
        }
        self.default_values = {
            "colors" : "black",
            "linestyles" : "-",
            "linewidths" : 2,
            "labels" : "_"
        }
        self.span_root_paths = {
            "labels" : True,
            "titles" : False,
            "linestyles" : True,
            "linewidths" : True,
            "colors" : True
        }

    def read_cmd_args(self):
        for argument in sys.argv[1:]:
            if "=" not in argument:
                self.attributes["root_paths"].append(argument)
            else:
                i=0
                while argument[i]!="=":
                    i+=1
                symbol_key = self.attributeSymbols[argument[:i]]
                if self.span_root_paths[symbol_key] == True:
                    self.attributes[symbol_key][self.attributes["root_paths"]
                        [-1]] = argument[i+1:]
                else:
                    self.attributes[symbol_key].append(argument[i+1:])

    def add_cmd_arg(self,arg_name,arg_symbol,span_root_paths=False):
        if arg_name[-1] != "s": #an "s" is added because we always assume there
                                #might be multiple arguments of this type
            arg_name+="s"
        if arg_symbol[0] != "-":
            arg_symbol = "-" + arg_symbol
        if span_root_paths:
            self.attributes[arg_name] = {}
            self.span_root_paths[arg_name] = True
        else:
            self.attributes[arg_name] = []
            self.span_root_paths[arg_name] = False
        self.attributeSymbols[arg_symbol] = arg_name

    def __str__(self):
        return repr(self.attributes)

    def get_attributes(self,name):
        if name[-1] != "s":
            name+="s"
        return self.attributes[name]

    def get_all_attributes(self):
        return self.attributes

    def get_attribute(self, type, key=0):
        if type[-1] != "s":
            type+="s"
        if self.span_root_paths[type]:
            attribute = self.attributes[type][key] if key in self.attributes[type] else self.default_values[type]
        else:
            attribute = self.attributes[type][key] if len(self.attributes[type])>0 else None
        return attribute

    def is_legend(self):
        legends={"_"}
        for key in self.attributes["labels"]:
            legends.add(self.attributes["labels"][key])
        if legends == {"_"}:
            return False
        else:
            return True
