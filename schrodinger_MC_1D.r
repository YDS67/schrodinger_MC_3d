# Решение уравнение Шредингера
# Задача Штурма-Лиувилля (методом конечных разностей) и метод Монте-Карло

h2m0 <- 0.0762 # eV*nm^2
me <- 0.067 #GaAs
a <- 4 #nm

U <- function(x){
	U0*(1 - exp(-(x-x1)^2/a^2) - exp(-(x-x2)^2/a^2))
}

L <- 60 #nm
x1 <- L/3
x2 <- 2*L/3
U0 <- 0.25 #eV
dx <- 0.1 #nm
Nx <- floor(L/dx)

E0 <- h2m0/me/dx^2
xn <- (1:Nx)*dx
Un <- U(xn)
un <- Un/E0

D <- h2m0/2/me
# D = dx^2/2/dtau, for 1D
dtau <- dx^2/2/D # 1/eV
p_evap <- Un*dtau
# p_evap >= 0 && p_evap <= 1

final_time <- 50 #1/eV
N_particles <- 2000
N_time_steps <- floor(final_time/dtau)
N_values_check <- 200
check_interval <- floor(N_time_steps / N_values_check)

filename <- paste("schr_FD_MC-", N_particles, "-", N_time_steps, "-", N_values_check, "-", Nx, sep="")

H <- diag(un+1)

H[row(H) - col(H) == 1] <- rep(-1/2, Nx-1)
H[row(H) - col(H) == -1] <- rep(-1/2, Nx-1)

# Находим собственные значения и собственные вектора матрицы Гамильтона H
eg <- eigen(H)

values <- eg$values*E0
vectors <- eg$vectors

#====================
# Monte-Carlo
#====================

energy_values <- NULL
tau_values <- NULL
evaporated_particles <- 0
t_current <- 0

particle_array <- rep(1,N_particles)
# Initial condition random
for(j in 1:N_particles){
	particle_array[j] <- sample(1:Nx,1)
}

#====================
# Цикл по времени - ЭТОТ ЦИКЛ НУЖНО УБРАТЬ В С++
#====================
for (jt in 1:N_time_steps) {
  # Цикл по частицам
  for (jn in 1:N_particles) {
    # Случайный скачок направо или налево
    particle_array[jn] <- particle_array[jn] + sample(c(-1,1),1)
	# Испарение с вероятностью p_evap
	if (particle_array[jn] >= 1 && particle_array[jn] <= Nx) {
	if (runif(1) < p_evap[particle_array[jn]]) {
		particle_array <- particle_array[-jn]
		particle_array <- c(particle_array, sample(particle_array, 1))
		evaporated_particles <- evaporated_particles+1
	}
	}
	else {
		particle_array <- particle_array[-jn]
		particle_array <- c(particle_array, sample(particle_array, 1))
	}
  }
  
  if (jt %% check_interval == 0) {
	tau_values <- c(tau_values, jt*dtau)
	# Вычисление энергии из скорости испарения
	if(N_particles-evaporated_particles > 0){
		energy <- log(N_particles/(N_particles-evaporated_particles))/(jt*dtau-t_current)
	}
	t_current <- jt*dtau
	evaporated_particles <- 0
	energy_values <- c(energy_values, energy)
  }
}

# Вводим функцию распеределения частиц
	particle_distribution <- rep(0, Nx)
	# Цикл по частицам
	for (jn in 1:N_particles) {
		if (particle_array[jn] >= 1 && particle_array[jn] <= Nx) {
			particle_distribution[particle_array[jn]] <- particle_distribution[particle_array[jn]] + 1
		}
	}
	particle_distribution <- particle_distribution/sqrt(sum(particle_distribution^2))

#==========================
# PLOTTING START
#==========================

nms <- c("U(x)",expression(paste(E[1]+abs(psi[1])^2, ", FD")), expression(paste(E[1]+abs(psi[1])^2, ", MC")))

x_max <- max(xn)
x_min <- min(xn)
x_length <- x_max-x_min
z_max <- 1.2*U0
z_min <- 0

number_of_lines <- 3

cols <- c("black", "#0092cc", "#ff3333")
ltys <- c(1, 1, 2)
pchs <- c(NA, NA, 16)

col_main <- "#282828"
lwd_main <- 5
axis_label_size <- 1.8
legend_font_size <- 1.5
tick_font_size <- 1.3
title_font_size <- 1.5

png(paste(filename, ".png", sep = ""), width = 3000, height = 2000, units = "px", pointsize = 45)
# par(mar = c(bottom, left, top, right))
par(mar = c(4.0, 4.2, 0.6, 0.6), mgp = c(2.5, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)

line_width <- 9
psz <- 1

plot(0, 0, pch = NA, # log = "y",
    xlim = c(x_min, x_max), ylim = c(z_min, z_max), 
    xlab = "x, nm", 
    ylab = "Energy, eV",
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

# grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)


lines(xn, Un, lwd = line_width, col = cols[1], lty = ltys[1])
abline(h=values[Nx], lwd = line_width, col = cols[2], lty = ltys[2])
lines(xn, values[Nx]+vectors[,Nx]^2/dx, lwd = line_width, col = cols[2], lty = ltys[2])

legend("topright", inset = 0.02, 
    legend = nms, 
    lwd = rep(line_width, number_of_lines), 
    lty = ltys,
    pch = pchs, 
    col = cols, 
    cex = legend_font_size,
    box.col = col_main)

N_values <- length(energy_values)
energy_correct <- mean(energy_values[floor(N_values/2):N_values])

#===================
# Plot the MC result
#===================

abline(h=energy_correct, lwd = line_width, col = cols[3], lty = ltys[3])
points(xn, energy_correct+particle_distribution^2/dx, col=cols[3], pch = 16)

dev.off()

#=====================
# Energy vs time plot
#=====================
nms <- c(expression(paste(E[1], ", FD")), expression(paste(E[1](tau), ", MC from evaporation rate")), expression(paste("<", E[1](tau), ">")))


x_max <- N_time_steps*dtau
x_min <- 0
x_length <- x_max-x_min
z_max <- U0
z_min <- 0

number_of_lines <- 3

cols <- c("#0092cc", "#ff3333", "#ff3333")
ltys <- c(1, 1, 2)
pchs <- c(NA, 16, NA)

col_main <- "#282828"
lwd_main <- 5
axis_label_size <- 1.8
legend_font_size <- 1.5
tick_font_size <- 1.3
title_font_size <- 1.5

png(paste(filename, "_energy.png", sep = ""), width = 3000, height = 2000, units = "px", pointsize = 45)
# par(mar = c(bottom, left, top, right))
par(mar = c(4.0, 4.2, 0.6, 0.6), mgp = c(2.5, 0.6, 0), lwd = lwd_main,
    bg = "#ffffff", col = col_main, col.lab = col_main, col.axis = col_main, col.main = col_main)

line_width <- 9
psz <- 1

plot(0, 0, pch = NA, # log = "y",
    xlim = c(x_min, x_max), ylim = c(z_min, z_max), 
    xlab = expression(paste(tau, ", ", eV^{-1})), 
    ylab = "Energy, eV",
    cex.axis = tick_font_size,
    cex.main = title_font_size,
    cex.lab = axis_label_size)

axis(1, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)
axis(2, col.ticks = col_main, lwd = lwd_main, lwd.ticks = lwd_main, col = col_main, labels = FALSE)

# grid(nx = 5, ny = 5, col = "lightgray", lty = 2, lwd = 0.5*line_width, equilogs = FALSE)

N_values <- length(energy_values)

abline(v=0.5*final_time, lwd = line_width, col = "darkgray", lty = 3)
abline(h=values[Nx], lwd = line_width, col = cols[1], lty = ltys[1])
points(tau_values, energy_values, col = cols[2], pch = pchs[2], cex = 1.5)
abline(h=energy_correct, lwd = line_width, col = cols[3], lty = ltys[3])

legend("topright", inset = 0.02, 
    legend = nms, 
    lwd = rep(line_width, number_of_lines), 
    lty = ltys,
    pch = pchs, 
    col = cols, 
    cex = legend_font_size,
    box.col = col_main)

dev.off()

write.table(data.frame(xn, vectors[,Nx], particle_distribution), paste(filename, "_density.dat", sep = ""), row.names = FALSE, col.names = FALSE)
write.table(data.frame(tau_values, energy_values, energy_values-values[Nx]), paste(filename, "_energy.dat", sep = ""), row.names = FALSE, col.names = FALSE)
