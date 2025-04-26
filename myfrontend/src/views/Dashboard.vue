<template>
  <div class="dashboard-container">
    <!-- 主要内容区域 -->
    <el-main>
      <el-row :gutter="20">
        <el-col :span="18">
          <div class="search-box">
            <el-input
              v-model="searchQuery"
              placeholder="请输入SMILES分子式"
              class="search-input"
              :disabled="!isLoggedIn"
            >
              <template #append>
                <el-button @click="handleSearch" :disabled="!isLoggedIn">
                  <el-icon><Search /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>

          <el-table
            v-loading="loading"
            :data="currentPageData"
            style="width: 100%"
            border
            stripe
            highlight-current-row
          >
            <el-table-column 
              type="index" 
              label="序号" 
              width="80" 
              align="center"
            />
            <el-table-column 
              prop="smiles" 
              label="SMILES分子式" 
              min-width="300"
            >
              <template #default="scope">
                <div 
                  class="clickable-cell"
                  @click="handleRowClick(scope.row)"
                >
                  {{ scope.row.smiles }}
                </div>
              </template>
            </el-table-column>
          </el-table>
          
          <!-- 添加分页组件 -->
          <div class="pagination-container">
            <el-pagination
              v-model:current-page="currentPage"
              v-model:page-size="pageSize"
              :page-sizes="[10, 20, 50, 100]"
              :total="smilesList.length"
              layout="total, sizes, prev, pager, next, jumper"
              @size-change="handleSizeChange"
              @current-change="handleCurrentChange"
            />
          </div>
        </el-col>
        
        <el-col :span="6" v-if="isLoggedIn">
          <el-card class="history-card">
            <template #header>
              <div class="card-header">
                <span>查询历史</span>
                <el-button type="text" @click="fetchQueryHistory">刷新</el-button>
              </div>
            </template>
            <el-timeline>
              <el-timeline-item
                v-for="item in queryHistory"
                :key="item.created_at"
                :timestamp="item.created_at"
                placement="top"
              >
                <div class="history-item" @click="searchQuery = item.smiles; handleSearch()">
                  {{ item.smiles }}
                </div>
              </el-timeline-item>
            </el-timeline>
          </el-card>
        </el-col>
      </el-row>
      
      <!-- 修改分子图片对话框 -->
      <el-dialog
        v-model="dialogVisible"
        title="分子结构图"
        :width="'70%'"
        center
      >
        <div class="image-container">
          <img :src="currentMoleculeImage" alt="分子结构图" v-if="currentMoleculeImage" class="molecule-image" />
        </div>
      </el-dialog>
    </el-main>
  </div>
</template>

<script setup>
import { ref, onMounted, watch } from 'vue'
import { useRouter } from 'vue-router'
import { ElMessage } from 'element-plus'
import { Search } from '@element-plus/icons-vue'
import axios from 'axios'

const router = useRouter()
const loading = ref(false)
const searchQuery = ref('')
const smilesList = ref([])
const dialogVisible = ref(false)
const currentMoleculeImage = ref('')
const isLoggedIn = ref(false)
const queryHistory = ref([])

// 添加分页相关的响应式变量
const currentPage = ref(1)
const pageSize = ref(10)
const currentPageData = ref([])

// 计算当前页的数据
const updateCurrentPageData = () => {
  const start = (currentPage.value - 1) * pageSize.value
  const end = start + pageSize.value
  currentPageData.value = smilesList.value.slice(start, end)
}

// 处理页码变化
const handleCurrentChange = (val) => {
  currentPage.value = val
  updateCurrentPageData()
}

// 处理每页条数变化
const handleSizeChange = (val) => {
  pageSize.value = val
  currentPage.value = 1
  updateCurrentPageData()
}

// 获取SMILES列表
const fetchSmilesList = async () => {
  try {
    loading.value = true
    const response = await axios.get('http://localhost:8000/api/smiles/')
    console.log('API Response:', response.data)
    
    if (response.data && response.data.prefixes) {
      smilesList.value = response.data.prefixes.map((smiles, index) => ({
        id: index + 1,
        smiles: smiles,
        created_at: new Date().toISOString()
      }))
      updateCurrentPageData()
    } else {
      smilesList.value = []
      currentPageData.value = []
      ElMessage.warning('数据格式不正确')
    }
  } catch (error) {
    console.error('Error fetching data:', error)
    ElMessage.error('获取数据失败：' + (error.response?.data?.detail || error.message))
    smilesList.value = []
    currentPageData.value = []
  } finally {
    loading.value = false
  }
}

// 处理搜索
const handleSearch = async () => {
  if (!isLoggedIn.value) {
    ElMessage.warning('请先登录后再进行搜索')
    return
  }
  
  if (!searchQuery.value) {
    ElMessage.warning('请输入要搜索的SMILES分子式')
    return
  }

  try {
    loading.value = true
    
    const token = localStorage.getItem('token')
    if (!token) {
      ElMessage.warning('登录已过期，请重新登录')
      router.push('/login')
      return
    }

    const response = await axios.get(`http://localhost:8000/api/query_smiles/?smiles=${encodeURIComponent(searchQuery.value)}`, {
      headers: {
        'Authorization': `Bearer ${token}`,
        'Accept': '*/*'
      },
      responseType: 'blob'
    })
    
    if (response.data) {
      const contentType = response.headers['content-type'] || 'image/png'
      const blob = new Blob([response.data], { type: contentType })
      const imageUrl = URL.createObjectURL(blob)
      currentMoleculeImage.value = imageUrl
      dialogVisible.value = true
      
      const unwatch = watch(dialogVisible, (newVal) => {
        if (!newVal) {
          URL.revokeObjectURL(imageUrl)
          unwatch()
        }
      })

      await fetchQueryHistory()
    }
  } catch (error) {
    console.error('Search error:', error)
    if (error.response?.status === 401) {
      localStorage.removeItem('token')
      isLoggedIn.value = false
      ElMessage.error('登录已过期，请重新登录')
      router.push('/login')
    } else {
      ElMessage.error('搜索失败：' + (error.response?.data?.detail || error.message))
    }
  } finally {
    loading.value = false
  }
}

// 处理行点击
const handleRowClick = async (row) => {
  console.log('Row clicked:', row)
  
  if (!isLoggedIn.value) {
    ElMessage.warning('请先登录后再查看分子结构图')
    return
  }
  
  try {
    loading.value = true
    const encodedSmiles = encodeURIComponent(row.smiles)
    
    const response = await axios.get(`http://localhost:8000/api/structure/?smiles=${encodedSmiles}`, {
      headers: {
        'Accept': '*/*'
      },
      responseType: 'blob'
    })
    
    if (response.data) {
      const contentType = response.headers['content-type'] || 'image/png'
      const blob = new Blob([response.data], { type: contentType })
      const imageUrl = URL.createObjectURL(blob)
      currentMoleculeImage.value = imageUrl
      dialogVisible.value = true
      
      const unwatch = watch(dialogVisible, (newVal) => {
        if (!newVal) {
          URL.revokeObjectURL(imageUrl)
          unwatch()
        }
      })
    }
  } catch (error) {
    console.error('Error fetching image:', error)
    ElMessage.error('获取分子结构图失败')
  } finally {
    loading.value = false
  }
}

// 获取查询历史
const fetchQueryHistory = async () => {
  if (!isLoggedIn.value) return
  
  try {
    const response = await axios.get('http://localhost:8000/api/query_history/')
    queryHistory.value = response.data.history
  } catch (error) {
    console.error('Error fetching history:', error)
    ElMessage.error('获取查询历史失败')
  }
}

// 检查登录状态
const checkLoginStatus = () => {
  const token = localStorage.getItem('token')
  isLoggedIn.value = !!token
  if (token) {
    axios.defaults.headers.common['Authorization'] = `Bearer ${token}`
  }
}

onMounted(() => {
  checkLoginStatus()
  fetchSmilesList()
  if (isLoggedIn.value) {
    fetchQueryHistory()
  }
})
</script>

<style scoped>
.dashboard-container {
  padding: 20px;
}

.search-box {
  margin-bottom: 20px;
}

.search-input {
  width: 100%;
}

.el-table {
  margin-top: 20px;
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

.clickable-cell {
  cursor: pointer;
  padding: 5px 0;
}

.clickable-cell:hover {
  color: #2c785c;
}

.history-card {
  margin-bottom: 20px;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.history-item {
  cursor: pointer;
  padding: 4px 0;
}

.history-item:hover {
  color: #2c785c;
}

.el-timeline {
  max-height: 600px;
  overflow-y: auto;
}

.pagination-container {
  margin-top: 20px;
  display: flex;
  justify-content: center;
}

.image-container {
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 20px;
}

.molecule-image {
  max-width: 100%;
  max-height: 70vh;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

:deep(.el-dialog__body) {
  padding: 20px;
  max-height: 80vh;
  overflow: auto;
}
</style> 