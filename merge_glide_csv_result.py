import os
import glob

base_path = "."  
output_file = "glide-dock_SP_chemdiv-stock1_01_400000.csv"
csv_files = []
for i in range(10):
    folder_name = f"glide-dock_SP_chemdiv-stock1_01_400000_split-10_{i}"
    csv_path = os.path.join(base_path, folder_name, f"glide-dock_SP_chemdiv-stock1_01_400000_split-10_{i}.csv")
    if os.path.exists(csv_path):
        csv_files.append(csv_path)

# 检查是否找到文件
if not csv_files:
    print("未找到CSV文件，请检查路径")
    exit()

# 按数字顺序排序文件（确保正确顺序）
csv_files.sort()
print(csv_files)
# 合并文件（跳过重复标题）
with open(output_file, 'w', newline='') as outfile:
    # 写入第一个文件的完整内容（含标题）
    with open(csv_files[0], 'r') as infile:
        outfile.write(infile.read())
    
    # 追加后续文件（跳过标题行）
    for file_path in csv_files[1:]:
        with open(file_path, 'r') as infile:
            # 跳过第一行（标题）
            next(infile)
            # 写入剩余内容
            outfile.write(infile.read())

print(f"合并完成！共合并 {len(csv_files)} 个文件到 {output_file}")
