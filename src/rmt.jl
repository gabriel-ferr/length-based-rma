function mean_recurrence_time(R::AbstractMatrix{Bool}; wmin = 1)
    N = size(R, 1)
    total = 0
    count = 0
    for j in axes(R, 2)
        last = 0                      # índice do último ponto recorrente na coluna (0 = nenhum ainda)
        for i in 1:N
            if R[i, j]
                if last > 0
                    w = i - last - 1  # comprimento da linha branca
                    if w ≥ wmin
                        total += w
                        count += 1
                    end
                end
                last = i
            end
        end
        # trecho após o último ponto recorrente é descartado (borda)
    end
    return count == 0 ? NaN : total / count
end