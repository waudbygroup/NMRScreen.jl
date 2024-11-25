"""
    @nested_observable outer inner
    
Creates an Observable that updates when either the outer Observable changes
or when the inner Observable (which is a field of outer's value) changes.
Automatically triggers initial setup by notifying the outer Observable.

Example:
    @nested_observable current_cocktail x
creates an Observable that updates when either current_cocktail changes
or when the current cocktail's x Observable changes.
"""
macro nested_observable(outer, inner)
    quote
        let
            # Get the initial value if possible
            try_initial = try
                $(esc(outer))[].$(esc(inner))[]
            catch
                Float64[]  # fallback initial value
            end
            @show try_initial
            # Create the result Observable
            result = Observable(try_initial)
            
            # Set up the nested listeners
            on($(esc(outer))) do outer_val
                on(getproperty(outer_val, $(QuoteNode(inner)))) do inner_val
                    result[] = inner_val
                end
                result[] = getproperty(outer_val, $(QuoteNode(inner)))[]
            end
            
            # Trigger initial setup
            notify($(esc(outer)))
            
            result
        end
    end
end


macro nested_observable_bool(outer, inner)
    quote
        let
            # Get the initial value if possible
            try_initial = try
                $(esc(outer))[].$(esc(inner))[]
            catch
                Bool[]  # fallback initial value
            end
            
            # Create the result Observable
            result = Observable(try_initial)
            
            # Set up the nested listeners
            on($(esc(outer))) do outer_val
                on(getproperty(outer_val, $(QuoteNode(inner)))) do inner_val
                    result[] = inner_val
                end
                result[] = getproperty(outer_val, $(QuoteNode(inner)))[]
            end
            
            # Trigger initial setup
            notify($(esc(outer)))
            
            result
        end
    end
end