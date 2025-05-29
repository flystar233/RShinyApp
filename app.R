library(shiny)
library(dplyr)
library(stringr)
library(ggplot2)

options(scipen=999,width = 200)
find_peaks <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector")
  }
  
  peaks <- integer(0)
  
  for (i in 2:(length(x) - 1)) {
    if (x[i] > x[i - 1] && x[i] > x[i + 1]) {
      peaks <- c(peaks, i)
    }
  }
  
  return(peaks)
}
kneedle <- function(x, y, decreasing, concave, sensitivity = 1) { # code from https://github.com/etam4260/kneedle/blob/main/R/kneedle.R
  if(length(x) == 0 || length(y) == 0) {
    stop("Make sure size of both inputs x and y are greater than 0")
  }
  if(typeof(x) == "list"|| typeof(y) == "list" || is.data.frame(x) ||
     is.data.frame(y) || is.array(x) || is.array(y) || is.matrix(x) || is.matrix(y)) {
    stop("Make sure both inputs x and y are vectors")
  }
  if(length(x) != length(y)) {
    stop("Make sure size of both inputs x and y are equal")
  }

  data <- matrix(unlist(list(x, y)), ncol = 2)
  data <- data[order(data[,1], decreasing = FALSE), ]

#  if(missing(decreasing)) {
#    if((data[(nrow(data)), 2] - data[1, 2]) >= 0) {
#      decreasing = FALSE
#    } else {
#      decreasing = TRUE
#    }
#  }
  if(missing(decreasing)) {
    fit <- lm(data[,2] ~ data[,1])
    if (coef(fit)[2] > 0) {
     decreasing = FALSE
     } else {
     decreasing = TRUE
     }
  }

  if(missing(concave)) {
    secondderiv <- diff(diff(data[, 2]) / diff(data[, 1]))
    if(mean(secondderiv) > 0) {
      concave = TRUE
    } else {
      concave = FALSE
    }
  }

  maxy <- max(y)
  miny <- min(y)
  maxx <- max(x)
  minx <- min(x)
  data[ ,1] <- (data[, 1]- min(data[, 1]))/(max(data[ ,1])- min(data[, 1]))
  data[ ,2] <- (data[, 2]- min(data[, 2]))/(max(data[ ,2])- min(data[, 2]))

  if(concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  } else if(!concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(!concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  }

  peak.indices <- find_peaks(differ)
  data <- cbind(data, differ)
  diffx = diff(data[, 1])
  T.lm.x.s <- sensitivity * mean(diffx)
  knee = NULL
  if(length(peak.indices) != 0){
    for(i in 1:length(peak.indices)) {
    T <- data[peak.indices[i] ,3] - (T.lm.x.s)

    for(j in (peak.indices[i]):if(i+1 < length(peak.indices)) peak.indices[i+1] else length(differ)) {
      if(differ[j] < T) {
        knee = peak.indices[i];
        break;
      }
    }
    if(!is.null(knee)) {
      break;
    }
    } 
  }

  x <- ((maxx - minx) * (data[knee, 1])) + minx
  y <- ((maxy - miny) * (data[knee, 2])) + miny

  return(c(as.numeric(x),as.numeric(y)))
}
ui <- fluidPage(
    titlePanel("Data Operation"),
    sidebarLayout(
        sidebarPanel(
            selectInput(
                inputId = "operation",
                label = "Choose an operation:",
                choices = c(
                    "交集(Intersect)", "并集(Union)",
                    "补集(Setdiff)", "去重(Unique)",
                    "排序(Sort)", "筛选(Filter by re)",
                    "筛选(Filter by sql)","云雨图(Raincloud plot)",
                    "时间戳转时间(Timestamp to Date)","相关性检验(Correlation)",
                    "拐点检测(knee Point plot)", "内容提换(Replace Content)"
                ),
                selected = "交集(intersect)"
            ),
            radioButtons("data", "Choose data:", choices = c("data1", "data2"), inline = TRUE),
            textAreaInput("data1", "data1", rows = 15, cols = 180),
            textAreaInput("data2", "data2", rows = 15, cols = 180),
            textInput("filter", "筛选条件（正则表达式）：", ".*"),
            textInput("filter2", "筛选条件(SQL):", "select * from df"),
            textAreaInput("replace_template", "内容提换模板(Template)", 
              value = '"{id}": {\n    "effective_date": "{date}"\n  },', 
              rows = 5, cols = 80),
        ),
        mainPanel(HTML('<textarea id="ta" class="form-control shiny-text-output"',
                       'style="resize:none;height:500px;" readonly></textarea>'),
                  br(),  
                  textOutput("text"),
                  verbatimTextOutput("COUNT"),
                  textOutput("text2"),
                  plotOutput("plot", width = "600px", height = "400px")
        )
    )
)

server <- function(input, output, session) {
    dataInput <- reactive({
        op <- input$operation
        data1 <- str_split(input$data1, "\n")[[1]]
        data2 <- str_split(input$data2, "\n")[[1]]
        filter <- input$filter
        data <- if (input$data == "data1") data1 else data2
        data_ <-if (input$data == "data1") data2 else data1
        result<-switch(
            op,
            "交集(Intersect)" = intersect(data1, data2),
            "并集(Union)" = union(data1, data2),
            "补集(Setdiff)" = setdiff(data, data_),
            "去重(Unique)" = unique(data),
            "排序(Sort)" = {
                if (all(str_detect(data, "^[0-9\\.]+$"))) {
                    sort(as.numeric(data))
                } else {
                    sort(data)
                }
            },
            "内容提换(Replace Content)" = {
            # data1 每行格式: id,date 或 id  date
            # 模板格式: "{id}": { "effective_date": "{date}" },
            template <- input$replace_template
            lines <- str_split(input$data1, "\n")[[1]]
            replaced <- sapply(lines, function(line) {
              # 自动判断分隔符：先尝试逗号，再尝试多个空格
              if (str_detect(line, ",")) {
                parts <- str_split(line, ",")[[1]]
              } else {
                # 用正则分割多个空格
                parts <- str_split(str_trim(line), "\\s{2,}|\\s+")[[1]]
              }
              if(length(parts) >= 2) {
                str_replace_all(template, c("\\{id\\}"=parts[1], "\\{date\\}"=parts[2]))
              } else {
                ""
              }
            })
            paste(replaced, collapse = "\n")
            },
            "时间戳转时间(Timestamp to Date)" = {
                timestamp <- as.numeric(data)/1000
                beijing_time <- as.POSIXct(timestamp, origin = "1970-01-01", tz = "Asia/Shanghai")
                new_york_time <- as.POSIXct(timestamp, origin = "1970-01-01", tz = "America/New_York")
                paste0("Beijing Time:\t", format(beijing_time, "%Y-%m-%d %H:%M:%S"), "\nNew York Time:\t", format(new_york_time, "%Y-%m-%d %H:%M:%S"))
            },
            "筛选(Filter by re)" = data[str_detect(data, filter)],
            "筛选(Filter by sql)" = {
                library(sqldf)
                #确定数据类型
                df_n10 <- read.table(text = data, header = TRUE, sep = "\t",nrows = 10)
                is_numeric <- sapply(df_n10, is.numeric)
                max_chars <- sapply(df_n10, function(x) max(nchar(as.character(x))))
                col_classes <- ifelse(max_chars > 15& is_numeric, "character", NA)
                
                df <- read.table(text=data,header = TRUE, sep = "\t",colClasses = col_classes)
                result<- sqldf(input$filter2)
                format_result<- paste(capture.output(print(result)), collapse = "\n")
            },
            "云雨图(Raincloud plot)" = {
                library(ggrain)
                dt <- read.table(text = data, header = FALSE, sep = "\t")
                ggplot(dt, aes(1,V1)) + 
                geom_rain(fill = "lightgray",alpha = .5,
                  point.args.pos = list(
                  position = position_jitter(width = 0.07)),
                  boxplot.args.pos = list(
                  width = 0.04, position = position_nudge(x = 0.13)),
                  violin.args.pos = list(
                  side = "r",
                  width = 0.7, position = position_nudge(x = 0.17)))+ 
                theme_classic()+xlab("")+ylab("")
            },
            "相关性检验(Correlation)" = {
                estimate_of_pearson<-""
                p_value_of_correlation_pearson<-""
                estimate_of_spearman<-""
                p_value_of_correlation_spearman<-""

                if (length(data1) == length(data2)) {
                    t_test <- t.test(as.numeric(data1),as.numeric(data2),paired=TRUE)
                    p_value_of_t_test <- format(t_test$p.value,scientific = TRUE)

                    wilcox <- wilcox.test(as.numeric(data1),as.numeric(data2),exact = F,paired=TRUE)
                    p_value_of_wilcox <- format(wilcox$p.value,scientific = TRUE)

                    correlation_pearson <- cor.test(as.numeric(data1), as.numeric(data2),method = 'pearson')
                    estimate_of_pearson <- correlation_pearson$estimate
                    p_value_of_correlation_pearson <-format(correlation_pearson$p.value,scientific = TRUE)

                    correlation_spearman<-cor.test(as.numeric(data1), as.numeric(data2),method = 'spearman')
                    estimate_of_spearman <- correlation_spearman$estimate
                    p_value_of_correlation_spearman <-format(correlation_spearman$p.value,scientific = TRUE)

            } else {
                    t_test <- t.test(as.numeric(data1),as.numeric(data2),paired=FALSE)
                    p_value_of_t_test <- format(t_test$p.value,scientific = TRUE)
                    wilcox <- wilcox.test(as.numeric(data1),as.numeric(data2),exact = F,paired=FALSE)
                    p_value_of_wilcox <- format(wilcox$p.value,scientific = TRUE)
                }
                paste0("T 检验(p value):\t",p_value_of_t_test,
                       "\nWilcox 秩和检验(p value):\t",p_value_of_wilcox,
                       "\nPearson 相关性检验(相关系数,p value):\t",estimate_of_pearson,"\t",p_value_of_correlation_pearson,
                       "\nSpearman 相关性检验(相关系数,p value):\t",estimate_of_spearman,"\t",p_value_of_correlation_pearson
                       )
            },
            "拐点检测(knee Point plot)" = {
                source("kpitu.R")
                in_data <- cbind(as.numeric(data1), as.numeric(data2))
                knee_point <- main_function(in_data,1)
                df <- data.frame(data1 = as.numeric(data1), data2 = as.numeric(data2))
                if(ncol(as.matrix(knee_point))==1)
                  {plot_knee<-ggplot() +
                   geom_point(data = df, aes(x = df[,1], y = df[,2])) +
                   geom_line(data = df, aes(x = df[,1], y = df[,2])) +
                   geom_point(aes(x = knee_point[1], y = knee_point[2],color = "red"), size = 4) +
                   annotate("text", x =knee_point[1] , y = knee_point[2], label = paste0("(",knee_point[1],",",knee_point[2],")"),
                            vjust=-1)+
                   labs(x="data1",y='data2',colour = "Knee point") +
                   guides(size = FALSE)
                   }else{
                  knee_point<-as.data.frame(knee_point)
                  plot_knee<-ggplot() +
                  geom_point(data = df, aes(x = df[,1], y = df[,2])) +
                  geom_line(data = df, aes(x = df[,1], y = df[,2])) +
                  geom_point(data = knee_point, aes(x = knee_point[,1], y = knee_point[,2],color = "red"), size = 4) +
                  annotate("text", x =knee_point[,1] , y = knee_point[,2], label = paste0("(",knee_point[,1],",",knee_point[,2],")"),
                           vjust=-1)+
                  labs(x="data1",y='data2',colour = "Knee point") +
                  guides(size = FALSE)
                    }
                plot_knee
            }
        )
        if(typeof(result) != "list"){
            paste(result, sep = "", collapse = "\n")
        }else {
            result
        }
        
    })

    output$ta<-renderText({
        req(!grepl("plot", input$operation)) 
        dataInput()
    })
    output$text <- renderText({ 
        "Count: " 
    })
    output$COUNT<-renderText({ 
        req(!grepl("plot", input$operation))
        data_2<-dataInput()
        sum(nzchar(str_split(data_2,"\n")[[1]]))
    })
     output$text2 <- renderText({ 
        "Plot: " 
    })
    output$plot <- renderPlot({
        req(grepl("plot", input$operation))
        print(dataInput())
    })
}

shinyApp(ui, server)
