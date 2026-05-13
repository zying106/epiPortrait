library(ggplot2)

# 创建存放 logo 的标准目录
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# 1. 生成六边形轮廓坐标
angles <- seq(0, 2 * pi, length.out = 7) + pi/6
hex <- data.frame(x = cos(angles), y = sin(angles))

# 2. 生成表观遗传学峰 (体现 4D 形状偏移)
x_val <- seq(-0.85, 0.85, length.out = 500)

# Peak 1: 尖峰 (代表 WT 或 Steep-Peak)
y1 <- dnorm(x_val, mean = -0.15, sd = 0.15) * 0.22 - 0.45
# Peak 2: 宽域、发生偏移 (代表 KO 或 Broad-Domain)
y2 <- dnorm(x_val, mean = 0.2, sd = 0.4) * 0.2 - 0.45

base_y <- -0.45
df1 <- data.frame(x = c(x_val, rev(x_val)), y = c(y1, rep(base_y, 500)))
df2 <- data.frame(x = c(x_val, rev(x_val)), y = c(y2, rep(base_y, 500)))

# 3. 使用 ggplot2 精密绘制
p <- ggplot() +
  # 蓝色的宽峰 (Broad-Domain)
  geom_polygon(data = df2, aes(x, y), fill = "#377EB8", alpha = 0.75) + 
  # 红色的尖峰 (Steep-Peak)
  geom_polygon(data = df1, aes(x, y), fill = "#E41A1C", alpha = 0.85) + 
  # 基因组基线
  geom_segment(aes(x = -0.85, xend = 0.85, y = base_y, yend = base_y), color = "#2C3E50", size = 1.2) +
  # 深色六边形外框
  geom_polygon(data = hex, aes(x, y), fill = NA, color = "#2C3E50", size = 3.5) +
  # 文字: epi (深沉稳重)
  annotate("text", x = 0, y = 0.35, label = "epi", size = 17, fontface = "bold", color = "#2C3E50", family = "sans") +
  # 文字: Portrait (亮眼斜体，体现多维画像)
  annotate("text", x = 0, y = 0.05, label = "Portrait", size = 12, fontface = "bold.italic", color = "#E41A1C", family = "sans") +
  coord_fixed() +
  theme_void() +
  theme(
    # 强制纯透明背景
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# 4. 导出为 300 PPI 高清透明 PNG
ggsave("man/figures/logo.png", plot = p, width = 5, height = 5, bg = "transparent", dpi = 300)
