classdef Cache < handle
    properties
        key
        elements
        head
        isFull
        cSize
    end
    
    methods
        
        function obj = Cache(cSize)
            obj.elements = cell(1,cSize);
            obj.key = cell(1,cSize);
            obj.head = 0;
            obj.isFull = false;
            obj.cSize = cSize;
        end
        
        function add(obj,key,val)
            obj.head = obj.head + 1;
            if obj.head > obj.cSize
                obj.head = 1;
                obj.isFull = true;
            end
            obj.key{obj.head} = key;
            obj.elements{obj.head} = val;
        end
        
        function [W] = get(obj,key)
            if obj.head == 0
                W = [];
            else
                n = obj.head;
                if obj.isFull;  n = obj.cSize; end
                for i = 1:n
                    if isequal(obj.key{i},key)
                        W = obj.elements{i};
                        return;
                    end
                end
                W = [];
            end
        end
    end
end
                    